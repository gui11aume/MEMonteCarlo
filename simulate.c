#include <math.h>
#include <stdio.h>

#include "mt.h"
#include "rbinom.h"

size_t
rgeom
(
   double p
)
{
   return ceil(log(runifMT())/log(1-p));
}


int
one_if_read_contains_MEM_seed
(
   size_t k,
   size_t N,
   size_t g,
   double p,
   double u
)
{

   size_t i;     // Size of the error-free segment.
   size_t m = 0; // Number of strictly masking threads.

   int sz = k;

   while (sz >= 0) {

      // Get the position of next error.
      i = rgeom(p);

      if (i > sz) {
         // If the error-free segment is longer than the read, the size
         // of the segment is k and the potentially masking threads
         // are not masking.
         
         // If the size of the read is shorter than the
         // minimum seed size, there is no seed at all.
         if (sz < g) return 0;

         // Otherwise, there is a seed if there is no surviving
         // strictly masking thread.
         return rbinom(m, pow(1-u,sz)) == 0;
      }
      else {
         // The next error is before the end of the read.
         // Get the number of strictly masking threads that survive
         // i nucleotides -- each survives with probability (1-u)^i
         size_t ssurv = rbinom(m, pow(1-u,i));

         // Get the number of potentially masking threads that survive
         // i+1 nucleotides -- each survives with probability (1-u)^i*u/3
         size_t psurv = rbinom(N-m, pow(1-u,i)*u/3.0);

         // If the segment has size greater than the minimum seed length
         // and there are no surviving threads then there is a seed.
         if (i >= g && ssurv == 0 && psurv == 0) return 1;

         // Otherwise, get the new number of strictly masking threads (m).
         // All the surviving potentially masking threads become strictly
         // masking. The others have a probability of u/3 of becoming
         // strictly masking.
         m = psurv + rbinom(N-psurv, u/3.0);

         // The new state is D/m/. Decrement 'k' and iterate.
         sz -= i+1; // NB: stop if last nucleotide is an error.
      }

   }

   // The end of the read is reached, no seed was found.
   return 0;

}


int main(void) {
   size_t tot = 0;
   for (int i = 0 ; i < 1000000 ; i++) {
      tot += one_if_read_contains_MEM_seed(50, 10000, 17, .01, .05);
   }
   fprintf(stdout, "%f\n", 1.0 - tot / 1000000.0);
}
