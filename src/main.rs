#![feature(stmt_expr_attributes)]
extern crate silva ;
extern crate time;
use std::env ;
use silva::* ;
use time::PreciseTime;
const SIGNBIT : i32 = 1<<31;
fn main() {
	let start = PreciseTime::now() ;
	let args: Vec<String> = env::args().collect();
	let exponent : i32 = args[1].parse().unwrap();
//	let exponent = 11 ;
	 let m = 10.0_f64.powi(exponent) as u64 ;
	 println!("{}",m) ;
  let mut beta = 0.00087 ;
//    if exponent = 18 then  beta:= 0.0033; 
  if exponent <= 7 { beta = 0.001 ; }
let alpha = beta * (exponent as f64 *10.0_f64.ln()).powi(3); // alpha = power(m,1/6) eliminates p2
 println!("alpha = {}",alpha);
let n = (alpha * (m as f64).cbrt()+0.5).floor() as u32;//trunc(alpha * power(m, 1 / 3) + 0.5);
let z = (10.0_f64.powf(exponent as f64 * 2.0 / 3.0) / alpha).floor() as u32 + 1 ;  //trunc(power(10, exponent * 2 / 3) / alpha ) + 1;
  println!("z = {} ", z);
  println!("n = {}", n);
  if n > z {
    println!("adjust beta");
    return;
    }
 let ll = (n+1) >> 1 ;
  //writeln('ll = ',ll) ;
let mut p : Box<[u32]> = vec![0;ll as usize + 3].into_boxed_slice() ; //: Vec<u32> = Vec::with_capacity(ll as usize + 3) ;//  Setlength(p,succ(ll)+2);
let mut mu : Box<[i32]> = vec![0;ll as usize + 3].into_boxed_slice() ; //Vec<i32> = Vec::with_capacity(ll as usize + 3) ;//  SetLength(mu, succ(LL)+2);
let mut pi : Box<[u32]> = vec![0;ll as usize + 3].into_boxed_slice() ; //Vec<u32> = Vec::with_capacity(ll as usize + 3) ;   //  Setlength(pi,succ(ll)+2);
  //init_mu(ll);
  let mut pix = 0 ;
  init_arrays(ll as usize,&mut mu,&mut pix,&mut pi, &mut p) ;
//println!("{:?}",p) ;  //for i in p do write(i,' ') ;
  //writeln();
//  println!("{:?}",mu) ;
  //for j in mu do write(j,' ') ;
  //writeln() ;
  //for i in pi do write (i,' ') ;
  println!("pix = {} ",pix) ;
// setlength(p,pix+2);
 let a = pi[(n as usize + 1) >> 1] as usize;
 let astar = pi[(int_sqrt(n as usize) + 1) >> 1] as usize;
  println!(" a = {}, astar = {}, p[a] = {} ",a,astar,p[a as usize]);
let lc = ((n as f64).log2()).floor() as u8 ;    //let lc = minm(26,trunc(log2(n)));
  println!(" lc = {} ", lc);
let  interval_length : usize = (1 << lc) as usize ;
let  num_intervals : usize = (z / interval_length as u32) as usize + 1 ;
  println!("# of intervals = {} ",num_intervals);
let mut interval_boundaries : Box<[u32]> = vec![0;num_intervals as usize + 2].into_boxed_slice() ;  //setlength(interval_boundaries , num_intervals+2);
let mut counter : Box<[i32]> = vec![0;interval_length as usize + 1].into_boxed_slice() ;   //  setlength(counter,interval_length+1) ;
let mut m1 : Box<[u32]> = vec![n;astar + 1].into_boxed_slice() ;   //  Setlength(m1, astar+1);
let mut phi : Box<[u64]> = vec![0;a + 2].into_boxed_slice() ;   //  Setlength(phi, succ(a)+1);
let mut t : Box<[u32]> = vec![0;a].into_boxed_slice() ;   //  Setlength(t, a);
let mut tt : Box<[u8]> = vec![0;a].into_boxed_slice() ;   //  Setlength(tt, a);
let mut d2 : Box<[u32]> = vec![0;a].into_boxed_slice() ;   //  Setlength(d2, a);
let mut offset : Box<[u32]> = vec![0;a + 2].into_boxed_slice() ;   //  SetLength(offset,succ(a)+1);
let mut block : Box<[bool]> = vec![false;n as usize + 4].into_boxed_slice() ;   //  SetLength(block,n+4) ;
let mut switch : Box<[bool]> = vec![false;a + 2].into_boxed_slice() ;   //  Setlength(switch, succ(a)+1);
  //switch := TBits.Create();
  //switch.Grow(succ(a));
 // switch.clearall();
 let mut  count : i64 = a as i64 - 1;
  for i in 0..num_intervals { interval_boundaries[i as usize] = 1 + (i * interval_length) as u32; }
  interval_boundaries[num_intervals] = z;
  //ord_leaves;
 ord_leaves(n,&mut count,&mu,m);
  // spec_leaves_beqc;
let mut  phi2 : i64 = (a as i64* (a as i64 - 1)) >> 1;
let mut  u = 10.0_f64.powf(exponent as f64 / 2.0).floor() as u32; //int_sqrt(m as usize);
  if u % 2 == 0 { u -= 1;}
  println!("u = {} ",u);
  let mut v = a;
  let mut w = u + 1;
  let mut s2bprimes = 0 ;
//  for prime in 1..astar { m1[prime as usize] = n; }
 for prime in astar..(a - 1) {
    s2b_init(prime,p[prime + 1],m,&mut t,n,&pi,a,&mut d2,&mut count,&mut tt) ;
    s2b_structured(prime,0,&mut s2bprimes,&mut d2,m,&p,&mut tt,n,&mut switch,&interval_boundaries,&mut count,&counter,&pi);
  }
//  setlength(t,0);
let  s1b_subst_const = 2;
for prime in 0..s1b_subst_const+1 { s1b_subst(prime,&p,n,&mu,m,&mut count) ; }
println!("count at start of main loop = {} ",count) ;
let mut p2primes = 0 ;
let mut endofprimes : usize = 0 ;
/// start of main loop
  for interval in 0..num_intervals {
   cnt_init(&mut counter, interval_length);
   for prime in 1..s1b_subst_const+1 {
    intervalclear(prime,&mut offset,&mut counter,interval_length,p[prime],lc) ; 
    }
    for prime in (s1b_subst_const + 1)..astar {
    intervalclear(prime,&mut offset,&mut counter,interval_length,p[prime],lc);
//m1b := m1[prime] ;   writeln('m1b = ', m1b) ;
s1b(prime,interval,&mut m1,n,p[prime + 1],m,&interval_boundaries,&mu,&mut count,&phi,&counter) ;
  phi[prime]+=(counter[interval_length - 1] & !SIGNBIT) as u64;
//        println!("interval = {}, count = {}, phi2 = {}, phi[prime] = {} ",interval, count , phi2, phi[prime as usize] );
    }
//          println!("interval = {}, count = {}, phi2 = {}, phi[prime] = {} ",interval, count , phi2, phi[prime as usize] );
      for prime in astar..(a - 1) {
     intervalclear(prime,&mut offset,&mut counter,interval_length,p[prime],lc);
       if switch[prime] {
 s2b_structured(prime,interval,&mut s2bprimes,&mut d2,m,&p,&mut tt,n,&mut switch,&interval_boundaries,&mut count,&counter,&pi);
count += (s2bprimes as u64 * phi[prime]) as i64 ;
 phi[prime]+= (counter[interval_length - 1] & !SIGNBIT) as u64;
    }
     else
      {
      	  endofprimes=prime; break;
      }
}
//      println!("interval = {}, count = {}, phi2 = {} ",interval, count , phi2 );
        for prime in  (endofprimes+1)..(a+1) { intervalclear(prime,&mut offset,&mut counter,interval_length as usize,p[prime],lc) ; }
   p2(interval,&mut p2primes,&mut u,&mut v,n,&mut w,&mut block,&p,m,&interval_boundaries,&mut phi2,&counter,a)  ;
      phi2 += phi[a] as i64 * p2primes as i64;
// println!("interval = {} , count = {}, phi2 = {}, p2primes = {} ",interval,count,phi2,p2primes);
 phi[a]+=(counter[interval_length - 1] & !SIGNBIT) as u64;
  
  }
  println!("phi2 = {} ",phi2) ;
   println!("prime count = {} ",count - phi2) ; 
   let end = PreciseTime::now();
    println!("{:?}", start.to(end)) ;            
}
