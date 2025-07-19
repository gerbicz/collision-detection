// Very fast search for collisions on a (large) 64 bits array.
// Written by Robert Gerbicz (2025). The original part is fun_collision_search(), the (trivial) compare64 is used only on a small
// array and on a last (very slow) check to see the correctness of this code (to find the repeated keys using quicksort).

// my long compilation line: g++ -m64 -O2 -fomit-frame-pointer -m64 -mtune=corei7 -march=corei7 -mavx2 -o collision collision.c -lm

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#if defined(WIN32) || defined(_WIN64)
	#define WIN32_LEAN_AND_MEAN

	#include <windows.h>
#else
	#include <fcntl.h>
	#include <unistd.h>
	#include <sys/resource.h>
#endif

#define MIN(a,b) ((a) < (b)? (a) : (b))
#define MAX(a,b) ((a) > (b)? (a) : (b))

typedef uint8_t uint8;
typedef uint32_t uint32;
typedef uint64_t uint64;

/* #define HAVE_PROF */
#ifdef HAVE_PROF
#define SHOW_PROF __attribute__((noinline))
#else
#define SHOW_PROF /* nothing */
#endif

/*---------------------------------------------------------------------*/
static void *
aligned_malloc(size_t len, uint32 align) {

	void *ptr, *aligned_ptr;
	unsigned long addr;

	ptr = malloc(len+align);

	 /* offset to next ALIGN-byte boundary */

	addr = (unsigned long)ptr;				
	addr = align - (addr & (align-1));
	aligned_ptr = (void *)((uint8 *)ptr + addr);

	*( (void **)aligned_ptr - 1 ) = ptr;
	return aligned_ptr;
}

/*---------------------------------------------------------------------*/
static void
aligned_free(void *newptr) {

	void *ptr;

	if (newptr == NULL) 
		return;
	ptr = *( (void **)newptr - 1 );
	free(ptr);
}

/*------------------------------------------------------------------*/
double
get_cpu_time(void) {

#if defined(WIN32) || defined(_WIN64)
	FILETIME create_time = {0, 0};
	FILETIME exit_time = {0, 0};
	FILETIME kernel_time = {0, 0};
	FILETIME user_time = {0, 0};

	GetThreadTimes(GetCurrentThread(),
			&create_time,
			&exit_time,
			&kernel_time,
			&user_time);

	return ((uint64)user_time.dwHighDateTime << 32 | 
	               user_time.dwLowDateTime) / 10000000.0;
#else
	struct rusage r_usage;

	#if 0 /* use for linux 2.6.26+ */
	getrusage(RUSAGE_THREAD, &r_usage);
	#else
	getrusage(RUSAGE_SELF, &r_usage);
	#endif

	return ((uint64)r_usage.ru_utime.tv_sec * 1000000 +
	               r_usage.ru_utime.tv_usec) / 1000000.0;
#endif
}

/*------------------------------------------------------------------------*/
static uint32 
get_rand(uint32 *rand_seed, uint32 *rand_carry) {
   
	/* A multiply-with-carry generator by George Marsaglia.
	   The period is about 2^63. */

	#define RAND_MULT 2131995753

	uint64 temp;

	temp = (uint64)(*rand_seed) * 
		       (uint64)RAND_MULT + 
		       (uint64)(*rand_carry);
	*rand_seed = (uint32)temp;
	*rand_carry = (uint32)(temp >> 32);
	return (uint32)temp;
}

/*------------------------------------------------------------------------*/
#define SORT_CACHE_SIZE (3 << 16)

#define MAX_HASHTABLE_BITS 16

typedef struct {
	uint32 num_sort;
	uint32 num_sort_alloc;
	uint64 *sort_key1;
	uint64 *sort_key2;
	uint32 *sort_data1;
	uint32 *sort_data2;

	uint8 * sort_cache;
} cpu_thread_data_t;

cpu_thread_data_t *
cpu_thread_data_init()
{
	cpu_thread_data_t *t = (cpu_thread_data_t *)calloc(1,
					sizeof(cpu_thread_data_t));

	t->sort_cache = (uint8 *)aligned_malloc(MAX(SORT_CACHE_SIZE, 
				sizeof(uint32) << MAX_HASHTABLE_BITS), 64);
	return t;
}

void
cpu_thread_data_free(cpu_thread_data_t * t)
{
	aligned_free(t->sort_cache);
	free(t->sort_key1);
	free(t->sort_key2);
	free(t->sort_data1);
	free(t->sort_data2);
	free(t);
}

static void 
grow_sort(cpu_thread_data_t * td, uint32 new_size)
{
	td->num_sort_alloc = new_size;
	td->sort_key1 = (uint64 *)realloc(td->sort_key1, 
				new_size * sizeof(uint64));
	td->sort_key2 = (uint64 *)realloc(td->sort_key2, 
				new_size * sizeof(uint64));
	td->sort_data1 = (uint32 *)realloc(td->sort_data1, 
				new_size * sizeof(uint32));
	td->sort_data2 = (uint32 *)realloc(td->sort_data2, 
				new_size * sizeof(uint32));
}
/*------------------------------------------------------------------------*/

#include <math.h> // for log(), rarely used.
#include <assert.h>
int compare64(const void* a, const void* b){
     uint64 a64=*((uint64*)a);
     uint64 b64=*((uint64*)b);
     
     if(a64==b64)return 0;
     else if(a64<b64) return -1;
     else return 1;
}

uint32 fun_collision_search(cpu_thread_data_t *td, uint32 key_bits, bool allow_loss){
// return by the number of near collisions, that is number of keys (counted with multiplicity) that is repeated in td.
// if allow_loss=true then it does not find a small 6-12 percent of repeated keys, but it is much faster.
    
    uint64 *key1=td->sort_key1;
    uint32 n=td->num_sort;

#define bucket_size 256 // should be power of two
#define bs1 255 // =bucket_size-1
#define SH 8 // =log2(bucket_size)

#define num_buckets 256 // should be power of two
#define HV 255          // =num_buckets-1
#define log2_num_buckets 8 // =log2(num_buckets)

// it is possible to use bucket_size!=num_buckets (for given n it is possible that for the optimal choice these are different)
    
#define d_extra_bits 4.5 //we use rounding, so it is possible that changing by say 0.99 it gives still the same rounded value for ilgo2,
                         //and gets the same speed for the same n value.

    assert(key_bits>=6+log2_num_buckets);
    assert(key_bits<=32+log2_num_buckets);
    
    uint32 i,j,offset[num_buckets];    
    uint32 ns=num_buckets+1+(n/bucket_size);
    uint32 last_bucket=num_buckets-1;
    uint32* nextblock=(uint32*)malloc(ns*sizeof(uint32));
    uint32* buckets=(uint32*)malloc(ns*bucket_size*sizeof(uint32));
    uint32 lenb[num_buckets];
    
    // cache friendly buckets, see the idea and code at https://sweet.ua.pt/tos/software/prime_sieve.html
    for(i=0;i<num_buckets;i++){offset[i]=i<<SH;lenb[i]=0;}
    for(i=0;i<n;i++){
        uint64 hv=key1[i]&HV;
        buckets[offset[hv]++]=key1[i]>>log2_num_buckets;
        if((offset[hv]&bs1)==0){
            nextblock[(offset[hv]-1)>>SH]=(++last_bucket);
            offset[hv]=last_bucket<<SH;
            lenb[hv]+=bucket_size;
        }
    }
    nextblock[last_bucket]=0;
    
    uint32 maxsize=1;
    for(i=0;i<num_buckets;i++){
        lenb[i]+=(offset[i]&bs1);
        maxsize=MAX(maxsize,lenb[i]);
    }

    uint32 *arr=(uint32*)malloc(maxsize*sizeof(uint32));
    uint32 block,p;
    
    int ilog2=MAX(5,(int)((double)log(maxsize+1)/log(2)+d_extra_bits));
    if(log2_num_buckets+ilog2>key_bits)ilog2=key_bits-log2_num_buckets;
    
    uint32 *T,*T2,*T3,*U,*U2,*U3;
    uint32 tsize=1u<<(ilog2-5);
    T=(uint32*)malloc(tsize*sizeof(uint32));
    T2=(uint32*)malloc(tsize*sizeof(uint32));
    T3=(uint32*)malloc(tsize*sizeof(uint32));
    
    uint32 cap=16,csize=0;
    uint64* C=(uint64*)malloc(cap*sizeof(uint64));
    // in C collect the possible repeated keys (it is possible before qsort that it contains also non-repeated key!)
    
    for(i=0;i<num_buckets;i++){
        
    if(lenb[i]==0)continue;

    uint32 cnt=0;
    for(j=0;j<tsize;j++)T[j]=0;
    
    if(!allow_loss){
       for(j=0;j<tsize;j++)T2[j]=0;
    }
    
    block=i;
    p=bucket_size*i;
    int my_shift,my_shift2=log2_num_buckets;
    
    ilog2=MAX(5,(int)((double)log(lenb[i]+1)/log(2)+d_extra_bits));
    if(log2_num_buckets+ilog2>key_bits)ilog2=key_bits-log2_num_buckets;
    my_shift2=0;
    uint32 hash_value,hash_value2=(1u<<ilog2)-1;
    
    for(;;){
        // see below, using the buckets build the first hashtable to find the possible repeated keys for key%num_buckets=i
        uint32 en=MIN(offset[i],p+bucket_size);
        
        if(!allow_loss){
           for(;p<en;){
               uint32 hv=buckets[p]&hash_value2;
               arr[cnt++]=buckets[p++];
        
               if(T2[hv>>5]&(1u<<(hv&31))){
                  T[hv>>5]|=(1u<<(hv&31));
               }
               else{
                  T2[hv>>5]|=(1u<<(hv&31));
               }
           }
        }
        else{
            for(;p<en;){
               uint32 hv=buckets[p]&hash_value2;
               arr[cnt++]=buckets[p++];        
               T[hv>>5]^=(1u<<(hv&31));
            }
        }
        if(p==offset[i])break;
        
        block=nextblock[block];
        if(block==0)break;
        p=block*bucket_size;
    }
    
#define maxit 20
    uint32 nsize[maxit];
    for(int it=0;it<maxit;it++){
    // If it%2==0 ---> hash data is in T array
    // If it%2==1 ---> hash data is in T3 array
        
    // Using hash table find possible repeated keys, if f(x) is not repeated then x is also not repeated.
    // For a hash function use f(x)=floor(x/2^k) mod (2^l) for different k,l values.
    // Using a table of size 2^l we can get the possible repeated keys, in another table store the possible repeated hash values.
    // in another pass extract the possible repeated keys.
    // In each ping-pong round we can reduce the size of the array (possible number of repeated keys).
        
        my_shift=my_shift2;
        hash_value=hash_value2;

        ilog2=MAX(5,(int)((double)log(cnt)/log(2)+d_extra_bits));
        if(ilog2>key_bits-log2_num_buckets)
           ilog2=key_bits-log2_num_buckets;
        
        if(it%2==0){
            U=T;U2=T2;U3=T3;
        }
        else{
            U=T3;U2=T2;U3=T;
        }
        
        my_shift2=my_shift+6;
        if(my_shift2+ilog2>key_bits-log2_num_buckets)my_shift2=0;
        
        hash_value2=(1u<<ilog2)-1;
        uint32 my_size=1u<<(ilog2-5);
                
        for(j=0;j<my_size;j++){
            U2[j]=0;
        }
        for(j=0;j<my_size;j++){
            U3[j]=0;
        }

        uint32 cnt2=0;
        
        if(allow_loss==false || it>0){
        for(j=0;j<cnt;j++){
            uint32 hv=(arr[j]>>my_shift)&hash_value;
            if(U[hv>>5]&(1u<<(hv&31))){
               uint32 hv2=(arr[j]>>my_shift2)&hash_value2;
               if(U2[hv2>>5]&(1u<<(hv2&31))){
                  U3[hv2>>5]|=1u<<(hv2&31);
               }
               else{
                  U2[hv2>>5]|=1u<<(hv2&31);
               }
               arr[cnt2++]=arr[j];
            }
        }}
        else{
        for(j=0;j<cnt;j++){
            uint32 hv=(arr[j]>>my_shift)&hash_value;
            if((U[hv>>5]&(1u<<(hv&31)))==0){
               uint32 hv2=(arr[j]>>my_shift2)&hash_value2;
               if(U2[hv2>>5]&(1u<<(hv2&31))){
                  U3[hv2>>5]|=1u<<(hv2&31);
               }
               else{
                  U2[hv2>>5]|=1u<<(hv2&31);
               }
               arr[cnt2++]=arr[j];
            }
        }
        }
        
        cnt=cnt2;
        nsize[it]=cnt;
        if(cnt==0 || it==maxit-1 || (it>=3 && nsize[it-3]==nsize[it])){
            while(csize+cnt>cap){
                cap+=cap/2;
                C=(uint64*)realloc(C,cap*sizeof(uint64));
            }
            for(j=0;j<cnt;j++)
                C[csize++]=(((uint64)arr[j])<<log2_num_buckets)+i;
            
            break;
        }
    }
    }
    
    // find the repeated terms in C (and store in the beginning of C exactly once)
    qsort(C,csize,sizeof(uint64),compare64);
    uint32 csize2=0;
    for(i=1;i<csize;i++)
        if(C[i]==C[i-1]&&(i==csize-1||C[i]!=C[i+1]))
            C[csize2++]=C[i];
    csize=csize2;
    if(csize==0)return 0;
    
// use buckets to find the repeated keys of td. first sort out the keys that are repeated, in D store the start position of the next(!) bucket
// then search the td's keys in this sorted array. (this is the known standard method).
    ilog2=MAX(5,(int)((double)log(csize+1)/log(2)+d_extra_bits));
    
    int ssize=1u<<(ilog2-5);
    uint32 hash_value=(1u<<ilog2)-1;
    int dsize=(1u<<ilog2)+1;
    uint32 *S=(uint32*)malloc(ssize*sizeof(uint32));
    uint32 *D=(uint32*)malloc(dsize*sizeof(uint32));
    uint64 *X=(uint64*)malloc(csize*sizeof(uint64));

    for(i=0;i<dsize;i++)D[i]=0;
    for(i=0;i<ssize;i++)S[i]=0;
    
    for(i=0;i<csize;i++){
        uint32 hv=C[i]&hash_value;
        D[hv+1]++;
        S[hv>>5]|=(1u<<(hv&31));
    }
    for(i=1;i<dsize;i++)D[i]+=D[i-1];
    for(i=0;i<csize;i++){
        X[D[C[i]&hash_value]++]=C[i];
    }

    uint32 near_collision=0;
    for(i=0;i<n;i++){
        uint32 hv=key1[i]&hash_value;
        if(S[hv>>5]&(1u<<(hv&31))){
           uint32 en=D[hv];
           for(j=(hv==0?0:D[hv-1]);j<en;j++)
               if(key1[i]==X[j]){
                  near_collision++;//do something, this key value is surely repeated.
                                   //so save data,key
                  break;
               }
        }
    }
    
    free(nextblock);
    free(buckets);
    free(arr);
    free(T);
    free(T2);
    free(T3);
    free(C);
    free(S);
    free(D);
    free(X);
    
    return near_collision;
}

/*------------------------------------------------------------------------*/
int main(int argc, char ** argv)
{
	uint32 i, j;
	uint32 key_bits = strtoul(argv[1], NULL, 10);
	uint32 num_sort = strtoul(argv[2], NULL, 10);
    bool allow_loss = (bool) atoi(argv[3]);
	double seconds;
	uint32 seed1 = 0x11111;
	uint32 seed2 = 0x22222;
	uint32 num_roots = 30000;
	uint32 pshift = 16;
	uint32 mask = (1 << pshift) - 1;
	uint32 id = 0;

	cpu_thread_data_t * ctx = cpu_thread_data_init();

	grow_sort(ctx, num_sort);

	for (i = j = 0; i < num_sort; i++) {
		uint64 key = (uint64)get_rand(&seed1, &seed2) << 32 |
				get_rand(&seed1, &seed2);
        
        //uint64 r1=get_rand(&seed1, &seed2);
        //uint64 r2=get_rand(&seed1, &seed2);
        //uint64 r3=get_rand(&seed1, &seed2);
        //key=(r1<<32)|r3;
        
        ctx->sort_key1[i] = key >> (64 - key_bits);
		ctx->sort_data1[i] = (id << pshift) |
				(get_rand(&seed1, &seed2) & mask);
		if (++j == num_roots) {
			j = 0;
			id++;
		}
	}
	ctx->num_sort = num_sort;

    seconds = get_cpu_time();
    uint32 near_collision=fun_collision_search(ctx,key_bits,allow_loss);
    double deltatime=get_cpu_time() - seconds;
    printf("On %u x %u-bit keys found %u near-collisions in %lf seconds, allowed_loss: %s.\n", 
	num_sort, key_bits,near_collision,deltatime,allow_loss?"Yes":"No");
    
    uint32 hits=0,num_diff=0;
    uint64* a=(uint64*)malloc(num_sort*sizeof(uint64));
    for(i=0;i<num_sort;i++)
        a[i]=ctx->sort_key1[i];
    qsort(a,num_sort,sizeof(uint64),compare64);
    for(i=0;i<num_sort;i++){
        if((i>0 && a[i]==a[i-1]) || (i+1<num_sort && a[i]==a[i+1]))hits++;
        if(i==0 || a[i]!=a[i-1]) num_diff++;
    }
    printf("Slow built-in quicksort reports %u near collisions. Number of different keys=%u.\n\n",hits,num_diff);
    free(a);
    
	cpu_thread_data_free(ctx);
	return 0;
}
