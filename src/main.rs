#![feature(stmt_expr_attributes)]
extern crate silva ;
extern crate chrono;
//use std::env ;
use std::io;
use silva::* ;
use chrono::* ;

const SIGNBIT : i32 = 1<<31;
const SUBSTITUTE : usize = 2 ;

fn main() {
'foo : loop {
println!("Please enter an integer from 1 to 17. The program will count the exact number of primes below this power of 10: ");
let mut exponent = String::new()  ;
io::stdin().read_line(&mut exponent).ok().expect("Failed to read line");
let exponent : i32 = exponent.trim().parse().ok().expect("Please enter an integer") ;
match exponent {
	1...17 => (), 
	_ => {println!("Not an integer between 1 and 17"); return ;},
	}
let start: DateTime<Local> = Local::now();
println!("{:?}",start) ;
 //	let args: Vec<String> = env::args().collect();
//	let exponent : i32 = args[1].parse().unwrap();
let m = 10.0_f64.powi(exponent) as u64 ;
//	 println!("{}",m) ;
let mut beta = 0.00087 ;
//    if exponent = 18 then  beta:= 0.0033; 
if exponent <= 7 { beta = 0.001 ; }
let alpha = beta * (exponent as f64 *10.0_f64.ln()).powi(3); // alpha = power(m,1/6) eliminates p2
// println!("alpha = {}",alpha);
let n = (alpha * (m as f64).cbrt()+0.5).floor() as usize;//trunc(alpha * power(m, 1 / 3) + 0.5);
let z = (10.0_f64.powf(exponent as f64 * 2.0 / 3.0) / alpha).floor() as usize + 1 ;  //trunc(power(10, exponent * 2 / 3) / alpha ) + 1;
//  println!("z = {} ", z);
//  println!("n = {}", n);
if n > z {
    println!("adjust beta");
    return;
}
 let mut ll = (n+1) >> 1 ;
 if exponent <= 5 { ll = (m as usize - 1) >> 1; }
let mut p : Box<[usize]> = vec![0;ll + 2].into_boxed_slice() ; 
let mut mu : Box<[isize]> = vec![1;ll + 2].into_boxed_slice() ;
let mut pi : Box<[usize]> = vec![0;ll + 2].into_boxed_slice() ;
let mut pix = 0 ;
initialize_arrays(ll,&mut mu,&mut pix,&mut pi, &mut p) ;
if exponent <= 5 {  
 println!("prime count = {} ", pix) ; continue;}
//println!("{:?}",p) ;  //for i in p do write(i,' ') ;
//  println!("{:?}",mu) ;
//  println!("pix = {} ",pix) ;
// setlength(p,pix+2);
let a = pi[(n as usize + 1) >> 1];
let astar = pi[(int_sqrt(n as usize) + 1) >> 1];
//  println!(" a = {}, astar = {}, p[a] = {} ",a,astar,p[a as usize]);
let lc = ((n as f64).log2()).floor() as u8 ;    //let lc = minm(26,trunc(log2(n)));
//  println!(" lc = {} ", lc);
let  interval_length = (1 << lc) as usize ;
let  num_intervals = (z / interval_length ) + 1 ;
println!("# of intervals = {} ",num_intervals);
let mut interval_boundaries : Box<[usize]> = vec![1;num_intervals + 1].into_boxed_slice() ;  //setlength(interval_boundaries , num_intervals+2);
let mut counter : Box<[i32]> = vec![0;interval_length].into_boxed_slice() ;   //  setlength(counter,interval_length+1) ;
let mut m1 : Box<[usize]> = vec![n;astar ].into_boxed_slice() ;   //  Setlength(m1, astar+1);
let mut phi : Box<[u64]> = vec![0;a + 1].into_boxed_slice() ;   //  Setlength(phi, succ(a)+1);
let mut t : Box<[usize]> = vec![0;a-1].into_boxed_slice() ;   //  Setlength(t, a);
let mut tt : Box<[u8]> = vec![0;a-1].into_boxed_slice() ;   //  Setlength(tt, a);
let mut d2 : Box<[usize]> = vec![0;a-1].into_boxed_slice() ;   //  Setlength(d2, a);
let mut offset : Box<[usize]> = vec![0;a + 1].into_boxed_slice() ;   //  SetLength(offset,succ(a)+1);
let mut block : Box<[bool]> = vec![false;n + 3].into_boxed_slice() ;   //  SetLength(block,n+4) ;
let mut switch : Box<[bool]> = vec![false;a + 1].into_boxed_slice() ;   //  Setlength(switch, succ(a)+1);
let mut  count : i64 = a as i64 - 1;
for i in 1..num_intervals { interval_boundaries[i] = 1 + (i * interval_length); }
interval_boundaries[num_intervals] = z;
ordinary_leaves(n,&mut count,&mu,m);
let mut  phi2 : i64 = (a as i64* (a as i64 - 1)) >> 1;
let mut  u  = 10.0_f64.powf(exponent as f64 / 2.0).floor() as usize; //int_sqrt(m as usize);
if u % 2 == 0 { u -= 1;}
//  println!("u = {} ",u);
let mut v = a;
let mut w = u + 1;
let mut s2bprimes = 0 ;
for prime in astar..(a - 1) {
    special_leaves_type_2_initialize(prime,p[prime + 1],m,&mut t,n,&pi,a,&mut d2,&mut count) ;
    special_leaves_type_2(prime,0,&mut s2bprimes,&mut d2,m,&p,&mut tt,n,&mut switch,&interval_boundaries,&mut count,&counter,&pi);
  }
//  setlength(t,0);
//let  s1b_subst_const = 2;
for prime in 0..SUBSTITUTE+1 { special_leaves_type_1_substitute(prime,&p,n,&mu,m,&mut count) ; }
println!("count at start of main loop = {} ",count) ;
let mut p2primes = 0 ;
let mut endofprimes : usize = a-2 ;
/// start of main loop
for interval in 0..num_intervals {
cnt_init(&mut counter, interval_length);
for prime in 1..SUBSTITUTE+1 {
    interval_clear(prime,&mut offset,&mut counter,interval_length,p[prime],lc) ; 
    }
for prime in (SUBSTITUTE + 1)..astar {
    interval_clear(prime,&mut offset,&mut counter,interval_length,p[prime],lc);
	special_leaves_type_1(prime,interval,&mut m1,n,p[prime + 1],m,&interval_boundaries,&mu,&mut count,&phi,&counter) ;
	phi[prime]+=(counter[interval_length - 1] & !SIGNBIT) as u64;
//        println!("interval = {}, count = {}, phi2 = {}, phi[prime] = {} ",interval, count , phi2, phi[prime as usize] );
}
//          println!("interval = {}, count = {}, phi2 = {}, phi[prime] = {} ",interval, count , phi2, phi[prime as usize] );
for prime in astar..(a - 1) {
     interval_clear(prime,&mut offset,&mut counter,interval_length,p[prime],lc);
     if switch[prime] {
 		special_leaves_type_2(prime,interval,&mut s2bprimes,&mut d2,m,&p,&mut tt,n,&mut switch,&interval_boundaries,&mut count,&counter,&pi);
 		count += (s2bprimes as u64 * phi[prime]) as i64 ;
 		phi[prime]+= (counter[interval_length - 1] & !SIGNBIT) as u64;
     }
     else
      {
      	endofprimes=prime; break;
      }
}
//      println!("interval = {}, count = {}, phi2 = {} ",interval, count , phi2 );
for prime in  (endofprimes+1)..(a+1) { interval_clear(prime,&mut offset,&mut counter,interval_length,p[prime],lc) ; } 
p2(interval,&mut p2primes,&mut u,&mut v,n,&mut w,&mut block,&p,m,&interval_boundaries,&mut phi2,&counter,a)  ;
phi2 += phi[a] as i64 * p2primes as i64;
// println!("interval = {} , count = {}, phi2 = {}, p2primes = {} ",interval,count,phi2,p2primes);
phi[a]+=(counter[interval_length - 1] & !SIGNBIT) as u64;
}
println!("phi2 = {} ",phi2) ;
println!("prime count for 10 ^ {} = {} ",exponent,count - phi2) ; 
let end: DateTime<Local> = Local::now();
println!("{:?}",end - start) ;
continue 'foo ;
}     
}
