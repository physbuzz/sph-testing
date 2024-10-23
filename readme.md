
This started as being focused on SPH, but instead let's see how far we can push CPU-based particle-particle interactions. 
If I'm running large code finding nearest neighbors, why not also investigate things like sphere packings? 
If I'm doing SPH in 2 and 3 dimensions, why not 4? etc. etc. :)

The dense code that I'm writing is likely very slow, but I want something that I can benchmark against.

- [x] VectorND class and basic tests
- [ ] DenseParticleGrid class and basic tests
- [ ] Allow periodic boundary conditions
- [ ] Implement basic fast graphics
- [ ] Test and benchmark with a jamming transition in 2D
- [ ] Test N-dimensional jamming configurations - I think up to 8 dimensions might be doable.
- [ ] Build jamming-nd standalone & plot that cool order parameter that I learned from Krauth.
- [ ] Implement basic SPH & do velocity plotting. Note that other things like Van der Waals gas would be easy too.

Questions that arise:
- Could I write iterators for the DenseParticleGrid and is it worth it? It would be nice if this naturally extended to sparse representations too.
- Is a Bresenham-Hypersphere algorithm easy to write? (loop over all cells that intersect the hypersphere (V,r) for r much greater than the cell size)
- 
