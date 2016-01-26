# silva-exact-prime-counting
A rust implementation of the algorithm of Tomas Oliveira e Silva for exact prime counting

When run from the command-line with 'cargo run --release' on a computer with nightly rust installed, this program prompts for an integer between 1 and 17. It will then output the exact number of primes below that power of 10, that is π(x) where x is 10^input. To work with stable rust the four step_by()'s need to be worked around with while loops and the doc markups removed.

For checking the answers there are online tables accessible from the home pages of Thomas R. Nicely and Tomas Oliveira e Silva and a print table forms an appendix of Hans Riesel's "Prime Numbers and Computer Methods for Factorization" 2/e.

The algorithm starts from Legendre's speeding up of the sieve of Eratosthenes. His method (see function rphi in the library) counts primes without actually generating them. It uses the inclusion-exclusion principle. The method was improved later in the nineteenth century by Meissel and then with the advent of electronic computers by D.H. Lehmer. These improvements involved reducing the a in the ϕ(x,a) term, which is the one that slows down the algorithm. A breakthrough was made by Lagaarias, Miller and Odlyzko in 1985, who found a complicated way to prune the binary tree that arises from the Legendre method. Deleglise and Rivat in 1996 found a further improvement and Oliveira e Silva's article, which gives a lot of help to the programmer and is available online, was published in 2006.
