#ifndef MY_TIMER_H_INCLUDED
#define MY_TIMER_H_INCLUDED

#include <ctime>

#define init_timer time_t T_TIMER_START, T_TIMER_END;\
                  clock_t C_TIMER_START, C_TIMER_END;
#define timeit(thing) T_TIMER_START = time(0); \
                      C_TIMER_START = clock(); \
                      thing; \
                      T_TIMER_END = time(0); \
                      C_TIMER_END = clock(); \
                      printf("Time = %.6f", 1.0*(T_TIMER_END - T_TIMER_START));\
                      printf(", %.6f\n", 1.0*(C_TIMER_END - C_TIMER_START)/CLOCKS_PER_SEC);

#endif // MY_TIMER_H_INCLUDED
