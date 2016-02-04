# silva-exact-prime-counting
A **Rust** implementation of the algorithm of Tomás Oliveira e Silva for exact prime counting

When run from the command-line with 'cargo run --release' on a computer with Rust installed, this program prompts for an integer between 1 and 17. It will then output the exact number of primes below that power (e) of 10, that is _π(x)_ where _x_ is 10^e. Then the user is prompted for another integer. To break the loop and return to the command prompt, Ctrl+C. 

For checking the answers there are online tables accessible from the home pages of Thomas R. Nicely and Tomás Oliveira e Silva and a print table forms an appendix of Hans Riesel's _Prime Numbers and Computer Methods for Factorization_ 2/e.

The algorithm starts from Legendre's speeding up of the sieve of Eratosthenes. His method (see function rphi in the library) counts primes without actually generating them. It uses the inclusion-exclusion principle. The method was improved later in the nineteenth century by Meissel and then with the advent of electronic computers by D.H. Lehmer. These improvements involved reducing the _a_ in the _ϕ(x,a)_ term, which is the one that slows down the algorithm. A breakthrough was made by Lagarias, Miller and Odlyzko in 1985, who found a way to prune the binary tree that arises from the Legendre method. Deléglise and Rivat in 1996 found a further improvement and Oliveira e Silva's article, which gives a lot of help to the programmer and is available online, was published in 2006. All these methods for the exact counting of primes are often classed as "combinatorial". Another class is "analytic" methods.

The Rust code is a translation of a Pascal program and far from idiomatic, to put it mildly. It takes about 60% of the time of the Pascal program. For 10^17 on an Intel Core i3-2100 CPU using a Linux x64 operating system it takes about ten minutes. It is hoped that Rust's support for concurrency can be used to improve this.
