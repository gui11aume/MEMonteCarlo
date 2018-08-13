#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mt.h"
#include "rbinom.h"


#define MEM_YES 1
#define MEM_NO 0


static size_t G = 17;
static size_t K = 50;
static double P = .01;
static double U = .05;

static double XI[51] = {0};
static double XIC[51] = {0};
static double ETA[51] = {0};
static double ETAC[51] = {0};

size_t
rgeom
(void)
{
   return ceil( log(runifMT()) / log(1-P) );
}


size_t
rpos
(
   const size_t m,
   const size_t i
)
{
   if (m == 0) return 1;
   const double mth_root_of_unif = pow(runifMT(), 1.0/m);
   return ceil( log(1-XI[i]*mth_root_of_unif) / log(1-U) );
}

void
fill_MEM_pos
(
   const size_t   N,
         size_t * pos
)
{

   size_t i;     // Size of the error-free segment.
   size_t m = 0; // Number of hard masks.

   // Note: stop if last nucleotide is an error.
   for (int sz = K ; sz > 0 ; sz -= (i+1)) {

      // Get size of the error-free segment.
      i = rgeom() - 1;

      if (i >= sz) {
         // If the error-free segment is longer than the read, resize
         // the error-free segment to 'sz' and ignore soft masks.
         
         // If the size of the read is shorter than the minimum
         // seed size 'g', then there is no on-target MEM seed.
         int long_enough = sz >= G;
         // Otherwise, there is an on-target MEM seed if there is
         // no survivng hard mask.
         int unmasked = rbinom(m, XIC[sz]) == 0;
         if (unmasked && long_enough) {
            // The segment is unmasked and long enough: there is a global
            // MEM seed. We need to find its position.
            size_t failpos = rpos(m, sz);
            size_t from = failpos < G ? K-sz + G : K-sz + failpos;
            for (int j = from ; j <= K ; j++) pos[j]++;
         }
         return;
      }
      else {
         // The next error is before the end of the read.
         // Get the number of hard masks that survive i
         // nucleotides -- each survives with probability (1-u)^i.
         size_t hsurv = rbinom(m, XIC[i]);

         // Get the number of soft masks that survive i+1
         // nucleotides -- each survives with probability (1-u)^i*u/3.
         size_t ssurv = rbinom(N-m, ETAC[i]);

         // If the segment has size greater than the
         // minimum seed length and there are no surviving
         // threads then there is an on-target MEM seed.
         int unmasked = hsurv == 0;
         int long_enough = i >= G;
         if (unmasked && long_enough) {
            int all_soft_masks_failed = ssurv == 0;
            size_t failpos = rpos(m, i);
            size_t from = failpos < G ? K-sz + G : K-sz + failpos;
            if (all_soft_masks_failed) {
               for (int j = from ; j <= K ; j++) pos[j]++;
               return;
            }
            else {
               size_t to = K-sz+i;
               for (int j = from ; j <= to ; j++) pos[j]++;
            }
         }

         // Otherwise, get the number of strictly masking threads (m).
         // All the surviving soft masks become strictly masking. The
         // hard masks have a probability u/3 of matching the error.
         // The soft masks that did not survive have a probability
         // pp = u/3 * xi / eta of matching the error. 
         double pp = 1 - (1 - U/3.0) / ETA[i];
         // The new state is down/m.
         m = ssurv + rbinom(m, U/3.0) + rbinom(N-m-ssurv, pp);

      }

   }

}


int main(void) {
   size_t pos[51] = {0};
   for (int i = 0 ; i <= 50 ; i++) {
      // Complement of 'xi'.
      XI[i]   = 1 - pow(1-U,i);
      XIC[i]  = 1 - XI[i];
      ETA[i]  = 1 - pow(1-U,i) * U/3.0;
      ETAC[i] = 1 - ETA[i];
   }
   for (int i = 0 ; i < 10000000 ; i++) {
      fill_MEM_pos(5, pos);
   }
   for (int i = 1 ; i <= 50 ; i++) {
      fprintf(stdout, "i:%d %.12f\n", i, 1.0 - pos[i] / 10000000.0);
   }
}
