#ifndef __PDM_TIMER_H__
#define __PDM_TIMER_H__

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Structure de mesure des temps d'execution
 *----------------------------------------------------------------------------*/

typedef struct {

  double  t_cpu;                   /* Temps CPU cumule */
  double  t_elapsed;               /* Temps elapsed cumule */
#if defined (PDM_HAVE_GETRUSAGE)
  double  t_cpu_u;                 /* Temps CPU utilisateur */
  double  t_cpu_s;                 /* Temps CPU system */
  double  t_cpu_debut;
  double  t_cpu_u_debut;
  double  t_cpu_s_debut;
#else
  clock_t t_cpu_debut;             /* Marque de debut de mesure
                                      du temps CPU */
#endif
  struct timeval t_elaps_debut;    /* Marque de debut de mesure
                                      du temps elapsed */
  int     indic;                   /* Indique si une mesure d'une tranche
                                      est en cours */
} PDM_timer_t;

/*============================================================================
 * Interfaces des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation d'un objet timer
 *
 * return
 *   timer
 *
 *----------------------------------------------------------------------------*/

PDM_timer_t *PDM_timer_create(void);

/*----------------------------------------------------------------------------
 * Reinitialisation des compteurs de temps
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

void PDM_timer_init(PDM_timer_t *timer);

/*----------------------------------------------------------------------------
 * Reprend la mesure du temps ecoule
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

void PDM_timer_resume(PDM_timer_t *timer);

/*----------------------------------------------------------------------------
 * Suspend la mesure du temps ecoule
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

void PDM_timer_hang_on(PDM_timer_t *timer);

/*----------------------------------------------------------------------------
 * Retourne le temps CPU en secondes
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

double PDM_timer_cpu(PDM_timer_t *timer);

/*----------------------------------------------------------------------------
 * Retourne le temps CPU utilisateur en secondes
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

double PDM_timer_cpu_user(PDM_timer_t *timer);

/*----------------------------------------------------------------------------
 * Retourne le temps CPU systeme en secondes
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

double PDM_timer_cpu_sys(PDM_timer_t *timer);

/*----------------------------------------------------------------------------
 * Retourne le temps elaps en secondes
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

double PDM_timer_elapsed(PDM_timer_t *timer);

/*----------------------------------------------------------------------------
 * Destruction d'un objet timer
 *
 * parameters :
 *   timer            <-- Timer
 *
 *----------------------------------------------------------------------------*/

void PDM_timer_free(PDM_timer_t *timer);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FICHIER_SEQ_H__ */
