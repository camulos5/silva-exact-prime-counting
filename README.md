# silva-exact-prime-counting
A rust implementation of the algorithm of Tomas Oliveira e Silva for exact prime counting

When run from the command-line with 'cargo run --release' on a computer with nightly rust installed, this program prompts for an integer between 1 and 17. It will then output the exact number of primes below that power of 10, that is π(x) where x is 10^input. To work with stable rust the four step_by()'s need to be worked around with while loops and the doc markups removed.
