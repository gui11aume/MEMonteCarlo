#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mt.h"
#include "rbinom.h"


size_t G;
size_t K;
double P;
double U;


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
   return ceil( log(1-(1-tXI[i])*Um) / log(1-U) );
}


void
fill_MEM_pos
(
   const size_t g,
   const size_t k,
   const double p,
   const double * txi,
   const size_t N,
         size_t * pos
)
{

   size_t i;     // Size of the error-free segment.
   size_t m = 0; // Number of hard masks.

   // Note: stop if last nucleotide is an error.
   for (int sz = k ; sz > 0 ; sz -= (i+1)) {

      // Get size of the error-free segment.
      i = rgeom(p) - 1;

      if (i >= sz) {
         // If the error-free segment is longer than the read, resize
         // the error-free segment to 'sz' and ignore soft masks.
         
         // If the size of the read is shorter than the minimum
         // seed size 'g', then there is no on-target MEM seed.
         // Otherwise, there is an on-target MEM seed if there is
         // no survivng hard mask.
         if (sz < g || rbinom(m, tXI[sz]) > 0) return;

         // Global MEM seed. Get the position.
         size_t failpos = rpos(m, sz);
         size_t from = failpos < g ? k-sz+g : k-sz+failpos;
         for (int j = from ; j <= k ; j++) pos[j]++;
         return;
      }
      else {
         // The next error is before the end of the read.
         // Get the number of hard masks that survive i
         // nucleotides -- each survives with probability (1-u)^i.
         size_t hsurv = rbinom(m, tXI[i]);

         // Get the number of soft masks that survive i+1
         // nucleotides -- each survives with probability (1-u)^i*u/3.
         size_t ssurv = rbinom(N-m, tXI[i]*u/3.0);

         // If the segment has size greater than the
         // minimum seed length and there are no surviving
         // threads then there is an on-target MEM seed.
         if (i >= g && hsurv == 0) {
            // All the hard mask failed. If the reads were shorter
            // than k, there would be shared MEM seed already so we
            // need to count it. If the soft masks fail after the
            // error-free segment, this is a "local" MEM seed that
            // is not present is extensions of the read. If the soft
            // masks fail within the error-free segment, the MEM
            // seed is "global" in the sense that it is present in
            // every extension of the read.
            size_t failpos = rpos(m, i);
            size_t from = failpos < g ? k-sz+g : k-sz+failpos;
            // We have the position where hard masks failed. Now
            // see if the MEM seed is local or global.
            double prob = pow( (1-tXI[i])/(1-tXI[i]*u/3.0), N-m);
            if (ssurv > 0 || runifMT() > prob) {
               // Local.
               size_t to = k-sz+i;
               for (int j = from ; j <= to ; j++) pos[j]++;
            }
            else {
               // Global.
               for (int j = from ; j <= k ; j++) pos[j]++;
               return;
            }
         }

         // Otherwise, get the number of strictly masking threads (m).
         // All the surviving soft masks become strictly masking. The
         // hard masks have a probability u/3 of matching the error.
         // The soft masks that did not survive have a probability
         // pp = u/3 * xi / eta of matching the error. 
         double pp = 1 - (1 - u/3.0) / (1 - tXI[i]*u/3.0);
         // The new state is down/m.
         m = ssurv + rbinom(m, u/3.0) + rbinom(N-m-ssurv, pp);

      }

   }

}


int main(void) {
   size_t pos[51] = {0};
   double txi[51] = {0};
   const double u = 0.05;
   const size_t k = 50;
   tXI = malloc((k+1) * sizeof(double));
   for (int i = 0 ; i <= k ; i++) {
      tXI[i] = pow(1-u,i);
   }

   for (int i = 0 ; i < 10000000 ; i++) {
      fill_MEM_pos(5, pos);
   }
   for (int i = 1 ; i <= 50 ; i++) {
      fprintf(stdout, "i:%d %f\n", i, 1.0 - pos[i] / 10000000.0);
   }
}
