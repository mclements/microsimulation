/*  Program to test the random number streams file:    RngStream.c   */

#include <stdio.h>
#include "RngStream.h"

#define NS 1000000

int main (void)
{
   double x;
   int i, k;
   RngStream g1, g2, g3;
   RngStream gar[NS];

   /* Create 3 parallel streams */
   g1 = RngStream_CreateStream ("Poisson");
   g2 = RngStream_CreateStream ("Cantor");
   g3 = RngStream_CreateStream ("Laplace");

   /* Generate 35 random integers in [5, 10] with stream g1 */
   for (i = 0;  i < 35; i++)
      k = RngStream_RandInt (g1, 5, 10);

   /* Generate 100 random reals in (0, 1) with stream g3 */
   for (i = 0;  i < 100;  i++)
      x = RngStream_RandU01 (g3);

   /* Restart stream g3 in its initial state and generate the same 100 
      random reals as above */
   RngStream_ResetStartStream (g3);
   for (i = 0; i < 100; i++)
      x = RngStream_RandU01 (g3);

   /* Send stream g3 to its next substream and generate 5  
      random reals in (0, 1) with double precision */
   RngStream_ResetNextSubstream (g3);
   RngStream_IncreasedPrecis (g3, 1);
   for (i = 0; i < 5; i++)
      x = RngStream_RandU01 (g3);

   /* Generate 100000 antithetic random reals in (0, 1) with stream g2 */
   RngStream_SetAntithetic (g2, 1);
   for (i = 0; i < 100000; i++)
      x = RngStream_RandU01 (g2);

   /* Delete the 3 streams */
   RngStream_DeleteStream (&g3);
   RngStream_DeleteStream (&g2);
   RngStream_DeleteStream (&g1);

   /* Create NS = 1000000 parallel streams */
   for (i = 0; i < NS; i++)
      gar[i] = RngStream_CreateStream ("");

   /* Generate 1000 random real in (0, 1) with stream 55555 in gar */
   for  (i = 0; i < 1000; i++)
      x = RngStream_RandU01 (gar[55554]);

   for  (i = 0; i < NS; i++)
      RngStream_DeleteStream (&gar[i]);

   return 0;
}
