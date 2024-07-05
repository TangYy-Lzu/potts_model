#ifndef parameter_h
#define parameter_h

#define J -1 // 铁磁为正，反铁磁为负
#define Q 3  // 状态数

#define TN 80  // 是温度的个数
#define B 120  // 120个退火，120个样本自发对称破缺到基态进行计算
#define N 1000 // 每个退火过程采样1000次

#define L 54 // 尺寸，diced为3倍数，unionjack为2倍数，centered-diced为6的整数倍

#define lattice diced // square,diced,union_jack
#define tcrit_up 0.7
#define tcrit_down 0.4

#define deltat 0.1
#define deltat_crit 0.01 // 相变点附近精确计算

#define WRITE false

#endif

/*
 * square 2:1-0.3
 *        3:
 * diced 2:
 *       3:0.4-0.7
 * unionjack 2:
 *           3:
 *
 */