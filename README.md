# BuzzHash

A locality sensitive hashing algorithm based on Desgupta, *et. al.* [A neural algorithm for a fundamental computing problem](http://dx.doi.org/10.1101/180471) (preprint).

An implementation of the algorithm described by Desgupta, *et. al.*, in the above reference and recently published in Science magazine [1]. The algorithm is based on the olfactory system of *Drosophila melanogaster*, though it is not a simulation of that system. It is, instead, a biomimetic abstraction applicable to the computer science problem of similarity search, but might well serve as an abstract model of neural organization applicable to other organisms and neural subsystems. In the last regard, the authors mention three possibilities: mouse olfaction, rat cerebellum, and rat hippocampus including entorhinal cortex, dentate gyrus, and hilar cells.

The algorithm is simple. Its input is a real vector, **x**, which in the fly's case has 50 components representing aggregate output of 50 types of olfactory receptor. The first step is to form **x-mean(x)** *i,e,* to subtract the mean of **x** from each of its elements. This mimics the biological step known as divisive normalization. The "mean centered" input is then multiplied by a sparse 1/0 matrix, **A**, which expands the number of components in **x-mean(x)**. In the fly's case, the 50 components are expanded to 2000, each of the 2000 values being an aggregate of about 6 original. The non-zero components of **A** are randomly chosen and, of course, once chosen are fixed. The final step is zeroize all but the largest few components of **A(x-mean(x))**, leaving a sparse vector to act as a tag. In the fly's case the top 5% are retained.

This package provides two algorithms, `sprand_fd` to form the matrix, **A**, and `buzzhash(A, x, topN)` to apply the algorithm to **x**, retaining the `topN` maximum values.

### References:

[1] S. Dasgupta, C. F. Stevens, and S. Navlakha (2017). [A neural algorithm for a fundamental computing problem.](http://science.sciencemag.org/content/358/6364/793.full?ijkey=aX3uts9Y4xqPE&keytype=ref&siteid=sci) Science, 358, 6364:793-796.

[2] Charles F. Stevens [A statistical property of fly odor responses is conserved across odors](http://www.pnas.org/content/113/24/6737.full)

[3] Charles F. Stevens [What the fly's nose tells the fly's brain](http://www.pnas.org/content/112/30/9460.full)

[future] John Myles White [MNIST.jl](https://github.com/johnmyleswhite/MNIST.jl#mnistjl)