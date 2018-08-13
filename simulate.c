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
   const double Um = pow(runifMT(), 1.0/m);
   return ceil( log(1-(1-pow(1-U,i))*Um) / log(1-U) );
}


int
fill_MEM_pos
(
   const size_t N
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
         int too_short = sz < G;
         // Otherwise, there is an on-target MEM seed if there is
         // no survivng hard mask.
         int masked = rbinom(m, pow(1-U,sz)) > 0;
         return (masked || too_short) ? MEM_NO : MEM_YES;
      }
      else {
         // The next error is before the end of the read.
         // Get the number of hard masks that survive i
         // nucleotides -- each survives with probability (1-u)^i.
         size_t hsurv = rbinom(m, pow(1-U,i));

         // Get the number of soft masks that survive i+1
         // nucleotides -- each survives with probability (1-u)^i*u/3.
         size_t ssurv = rbinom(N-m, pow(1-U,i)*U/3.0);

         // If the segment has size greater than the
         // minimum seed length and there are no surviving
         // threads then there is an on-target MEM seed.
         int unmasked = hsurv == 0 && ssurv == 0;
         int long_enough = i >= G;
         if (unmasked && long_enough) return MEM_YES;

         // Otherwise, get the number of strictly masking threads (m).
         // All the surviving soft masks become strictly masking. The
         // hard masks have a probability u/3 of matching the error.
         // The soft masks that did not survive have a probability
         // pp = u/3 * xi / eta of matching the error. 
         double pp = 1 - (1 - U/3.0) / (1 - pow(1-U,i)*U/3.0);
         // The new state is down/m.
         m = ssurv + rbinom(m, U/3.0) + rbinom(N-m-ssurv, pp);

      }

   }

   return MEM_NO;

}


int main(void) {
   double tot = 0.0;
   for (int i = 0 ; i < 10000000 ; i++) {
      tot += fill_MEM_pos(5);
   }
   fprintf(stdout, "%.12f\n", 1.0 - tot / 10000000.0);
}
