/*
 *  Bit rank and select, taken from CQF (Pandey et al.)
 */

#ifndef RANKSELECT_H
#define RANKSELECT_H

int tzcnt(uint64_t val);
int popcnt(uint64_t val);
uint64_t bitrank(uint64_t val, uint64_t pos);
uint64_t bitselect(uint64_t val, uint64_t rank);

#endif //RANKSELECT_H
