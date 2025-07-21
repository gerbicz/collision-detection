Usage: ./collision k n b

to find (near) collisions on n number of keys, each key has k bits. If b=0 then we find all such near collisions, for b=1 it does not find a small 6-12 percent of repeated keys, but it is much faster. On my core-i7 a sample run:

gerbicz@gerbicz-MS-7972:~$ ./collision 39 20000000 0
On 20000000 x 39-bit keys found 1314 near-collisions in 0.170347 seconds, allowed_loss: No.
Slow built-in quicksort reports 1314 near collisions. Number of different keys=19999343.

gerbicz@gerbicz-MS-7972:~$ ./collision 39 20000000 1
On 20000000 x 39-bit keys found 1234 near-collisions in 0.162624 seconds, allowed_loss: Yes.
Slow built-in quicksort reports 1314 near collisions. Number of different keys=19999343.

gerbicz@gerbicz-MS-7972:~$
