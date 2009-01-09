#ifndef __SOLVE_AX_B_4_H__
#define __SOLVE_AX_B_4_H__

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*----------------------------------------------------------------------------
 * Solve Ax = B for a 4x4 system.
 *
 * parameters:
 *   a                  <-- matrix
 *   b                  <-- right hand side
 *   x                  --> solution of Ax = b
 *
 * returns: 0 if solution was found, 1 in case of error (zero pivot)
 *----------------------------------------------------------------------------*/

//int solve_ax_b_4(double a[4][4], double  *restrict b, double  *restrict x);
int solve_ax_b_4(double a[4][4], double  * b, double  * x);
#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __BAR_COORDS_H__ */
