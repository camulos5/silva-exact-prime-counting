//#![feature(stmt_expr_attributes)]
extern crate silva ;
extern crate chrono ;
extern crate bit_vec ;
extern crate itertools ;

use std::io;
use bit_vec::BitVec;
use silva::* ;
use chrono::* ;
use itertools::* ;

const SIGNBIT : i32 = 1<<31;
const SUBSTITUTE : usize = 2 ;

fn main() {
'foo : loop {
println!("Please enter an integer from 1 to 17. The program will count the exact number of primes below this power of 10: ");
let mut exponent = String::new()  ;
io::stdin().read_line(&mut exponent).ok().expect("Failed to read line");
let exponent : u32 = exponent.trim().parse().ok().expect("Please enter an integer") ;
match exponent {
	1...17 => (), 
	_ => {println!("Not a valid input: integer between 1 and 17"); return ;},
	}
let start: DateTime<Local> = Local::now();
println!("{:?}",start) ;
let m = 10u64.pow(exponent) ;
let mut beta = 0.00087 ;
//    if exponent = 18 then  beta:= 0.0033; also needs posssibly quite extensive re-typing and takes over an hour
if exponent <= 7 { beta = 0.001 ; }
let alpha = beta * (exponent as f64 *10.0_f64.ln()).powi(3); 
let n = (alpha * (m as f64).cbrt()+0.5).floor() as usize;
let z = (10.0_f64.powf(exponent as f64 * 2.0 / 3.0) / alpha).floor() as usize + 1 ;  
if n > z {
    println!("adjust beta");
    return;
}
let mut ll = (n+1) >> 1 ;
if exponent <= 5 { ll = (m as usize - 1) >> 1; }
let mut primes : Vec<usize> = vec!(1;ll+2) ;
let mut mu : Vec<isize> = vec![1;ll + 2];
let mut pi : Vec<usize> = vec![0;ll + 2];
let pix = initialize_arrays(ll,&mut mu,&mut pi, &mut primes) ;
if exponent <= 5 {  
 println!("prime count = {} ", pix) ;
let end: DateTime<Local> = Local::now();
println!("{:?}",end - start) ; 
 continue;
 }
primes.truncate(pix+1) ;
let a = pi[(n + 1) >> 1];
let astar = pi[(int_sqrt(n) + 1) >> 1];
let lc = ((n as f64).log2()).floor() as u8 ;   
let  interval_length = (1 << lc) as usize ;
let  num_intervals = (z / interval_length ) + 1 ;
let mut interval_boundaries : Vec<usize> = vec![1;num_intervals + 1];
let mut initial : Vec<i32> = vec![0;interval_length];
let mut m1 : Vec<usize> = vec![n;astar ];
let mut phi : Vec<u64> = vec![0;a + 1];
let mut t : Vec<usize> = vec![0;a-1];
let mut tt : Vec<u8> = vec![0;a-1];
let mut d2 : Vec<usize> = vec![0;a-1];
let mut offsets : Vec<usize> = vec![0;a + 1];
let mut block : BitVec = BitVec::from_elem(n+3,false);
let mut switch : Vec<bool> = vec![false;a + 1];
let mut  count : i64 = a as i64 - 1;
let intervals = (0..num_intervals).collect::<Vec<usize>>() ;
for i in 1..num_intervals { interval_boundaries[i] = 1 + (i * interval_length); }
interval_boundaries[num_intervals] = z;
ordinary_leaves(n,&mut count,&mu,m);
let mut  phi2 : i64 = (a as i64* (a as i64 - 1)) >> 1;
let mut u = match exponent % 2 {
 0 =>  10usize.pow(exponent/2) - 1,
 _  => 10.0_f64.powf(exponent as f64 / 2.0).floor() as usize,
 }; //int_sqrt(m as usize);
if u % 2 == 0 { u -= 1;}
let mut v = a;
let mut w = u + 1;
let mut p2primes = 0 ;
let mut s2bprimes = 0 ;
let mut endofprimes : usize = a-2 ;
 for index in 0..(SUBSTITUTE+1) {
 	count -= special_leaves_type_1_substitute(index,&primes,n,&mu,m) ; 
  }  ;
for index in astar..(a - 1) {
    special_leaves_type_2_initialize(index,primes[index + 1],m,&mut t,n,&pi,a,&mut d2,&mut count) ;
    special_leaves_type_2(index,0,&mut s2bprimes,&mut d2,m,&primes,&mut tt,n,&mut switch,&interval_boundaries,&mut count,&initial,&pi);
  }
	initial.iter_mut().into_rc().enumerate().map( |(i,e)| {*e = (i as i32 +1) & !(i as i32) } ).collect_vec() ;

// start of main loop
for interval in intervals {
let mut	 counter = &mut initial.clone() ;
for index in 1..SUBSTITUTE+1 {
  interval_clear(index,&mut offsets,&mut counter,interval_length,primes[index]) ; 
    }
for index in SUBSTITUTE+1..astar {
    interval_clear(index,&mut offsets,&mut counter,interval_length,primes[index]) ;
	special_leaves_type_1(index,interval,&mut m1,n,primes[index + 1],m,&interval_boundaries,&mu,&mut count,&phi,&counter) ;
	phi[index] += (counter[interval_length - 1] & !SIGNBIT) as u64;
}
for index in astar..a-2 {
    interval_clear(index,&mut offsets,&mut counter,interval_length,primes[index]) ;
     if switch[index] {
 		special_leaves_type_2(index,interval,&mut s2bprimes,&mut d2,m,&primes,&mut tt,n,&mut switch,&interval_boundaries,&mut count,&counter,&pi);
 		count += (s2bprimes as u64 * phi[index]) as i64 ;
 		phi[index] += (counter[interval_length - 1] & !SIGNBIT) as u64;
     }
     else
      {
      	endofprimes=index; break;
      }
}
for index in endofprimes+1..a+1 {
	interval_clear(index,&mut offsets,&mut counter,interval_length,primes[index]) ; } 
p2(interval,&mut p2primes,&mut u,&mut v,n,&mut w,&mut block,&primes,m,&interval_boundaries,&mut phi2,&counter,a)  ;
phi2 += phi[a] as i64 * p2primes as i64;
phi[a] += (counter[interval_length - 1] & !SIGNBIT) as u64;
}
//end of main loop
println!("prime count for 10 ^ {} = {} ",exponent,count - phi2) ; 
let end: DateTime<Local> = Local::now();
println!("{:?}",end - start) ;
continue 'foo ;
}     
}
