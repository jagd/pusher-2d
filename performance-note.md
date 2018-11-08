#    Performance Benchmark

##   Intel(R) Core(TM) i7-4700MQ CPU @ 2.40GHz
###  Win10 VC2017-Release x64
#### With Gamma Correction
LeapFrog pusher: 16.9177s
Boris pusher: 23.4777s
A-Phi pusher: 10.0256s
RK4 pusher: 130.429s
#### Without Gamma Correction
LeapFrog pusher: 13.9076s
Boris pusher: 23.4611s
A-Phi pusher: 7.89667s
RK4 pusher: 131.201s

## AMD A10 PRO-7800B R7 @ 3.5 GHz
###  Win10 VC2017-Release x86 with Gamma Correction
LeapFrog pusher: 17.0773s
Boris pusher: 20.706s
A-Phi pusher: 13.9485s
RK4 pusher: 156.18s

## Intel(R) Xeon(R) CPU E5-2690 0 @ 2.90GHz
## Ubuntu 18.04 gcc 7.3 x64 with Gamma Correction
LeapFrog pusher: 11.2416s
Boris pusher: 14.8186s
A-Phi pusher: 8.80014s
RK4 pusher: 90.7538s

## Ubuntu 18.04 icc 19 x64 with Gamma Correction
LeapFrog pusher: 11.0553s
Boris pusher: 14.209s
A-Phi pusher: 8.72603s
RK4 pusher: 85.2278s
