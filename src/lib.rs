//#![feature(step_trait)]

extern crate bit_vec ;
extern crate itertools ;
extern crate core ;

use bit_vec::BitVec;
use std::cmp ;
use itertools::Itertools ;
//use core::iter::Step ;

const SIGNBIT : i32 = 1<<31;

pub fn modneg( x : i64, y : usize) -> usize {
	(( x % y as i64) + y as i64) as usize % y
	}


pub fn int_sqrt(n :usize) -> usize {
	((n as f64).sqrt()).floor() as usize
}  

#[inline]
pub fn cnt_query(mut pos : usize, counter : &[i32]) -> u32 {
	let mut acc = counter[pos] & !SIGNBIT  ;
	pos += 1 ;
	pos &= pos - 1 ;
	while pos  != 0 {	
		acc += counter[pos - 1] & !SIGNBIT;
		pos &= pos - 1;
	}
acc as u32
}

#[inline]
pub fn cnt_update(mut pos : usize, counter : &mut[i32], interval_length : usize ) {
    counter[pos] |= SIGNBIT;
    loop {
      counter[pos] -= 1;
      pos |= pos + 1; 
      if pos >=interval_length {break;}
    }
    }

#[inline]
pub fn interval_clear(index : usize , offsets : &mut[usize],  counter : &mut[i32], interval_length : usize, prime : usize) {
//   (offsets[index]..interval_length).step(prime).foreach(|j| {    	if counter[j] > 0 { cnt_update( j, counter, interval_length ) ; } }) ; much slower
  
   let mut j  = offsets[index] ;
   while j < interval_length {
   	if counter[j] > 0 { cnt_update( j, counter, interval_length ) ; }
   	j += prime ;
   	}
   offsets[index] = modneg( offsets[index] as i64 - interval_length as i64 , prime );
   }


 pub fn initialize_arrays(ll : usize, mu : &mut[isize], pi : &mut [usize],  primes : &mut [usize]) -> usize {
 	 	for j in  2..mu.len() {
 		if mu[j] == 1 {
 			let mut i = j ; while i <= ll  {
 			 mu[i] = match mu[i] {
 				1 => 1 - 2 * j as isize,
 				_ => -mu[i] ,
 			} ;	
 		i += 2 * j - 1 ;	}
 		}
 	} 
 	let mut j = 2 ; while j * j <= ll << 1 {
 		if mu[j] == 1- 2*j as isize {
 			let mut i = 2 * j * j - 2 * j + 1 ; 
 			while i <= ll  {
 			mu[i] = 0 ;	
 			i += 4 * j * j - 4 * j + 1 ;
 			}
 		}
 		j += 1 ;
 	}
primes[1] = 2 ;
let mut pix = 1 ;
for (i , elem) in mu.iter().enumerate().dropping(2) {
if *elem == 1 - 2 * i as isize {
 	pix += 1;
 	primes[pix] = 2 * i as usize -1 ;
}
pi[i as usize]=pix;	
 }
//println!("{:?}",mu) ;
 pix
 } 
 
/* 
 pub fn ordinary_leaves(n : usize, count : &mut i64, mu : &[isize], m : u64) {
  for i in (1..n + 1).step(2) {
  	let term = (mu[(i+1) >> 1]).signum() as i64;
  	*count +=  term * (m as i64/ i as i64) ;
 }
for i in (2..n + 1).step(4) {
	let term = (mu[((i >> 1)  + 1) >> 1]).signum() as i64 ;
	*count -= term * (m as i64 / i as i64) ;
}
}       
 */
  pub fn ordinary_leaves(n : usize, count : &mut i64, mu : &[isize], m : u64) {
let mut it = (1..n+1).filter(|&i| i%4 != 0) ;
it.foreach(|i| {if i % 2 == 1 
		{let term = (mu[(i+1) >> 1]).signum() as i64;
  	*count +=  term * (m as i64/ i as i64) ;}
		 else { let term = (mu[((i >> 1)  + 1) >> 1]).signum() as i64 ;
	*count -= term * (m as i64 / i as i64) ;} })
}
 

 pub fn rphi(term: i64,  mut a : usize, bit : i8 , primes: &[usize],  total : &mut i64 )  {
	loop {
		if a == 1 {
			*total +=  bit as i64 * term ; 
			break;
			}
		else if term < primes[a] as i64 {
			*total += bit as i64 ;
			break;
				}
		else {
			a -= 1;
			rphi(term / (primes[a] as i64), a, -bit, primes, total) ;
					}
				}
		}
 
 
pub fn special_leaves_type_1_substitute(b : usize, primes : &[usize], n : usize, mu : &[isize], m : u64) -> i64 {
let pp = primes[b + 1]    ; let mut acc = 0 ;
let mut j = cmp::max((n / pp) , pp ) as usize +1 ;
if j % 2 == 0 { j += 1;} 
let mut i = j ; while i <= n {
//for i in (j..(n+1)).step_by(2) {
let muval1 = (mu[ ( i+1) >> 1 ] ).signum() ;
if (mu[(i + 1) >> 1]).abs() > pp as isize && muval1 != 0 {
let term = (m / (i * pp) as u64) as i64;
let mut total : i64 = 0;
rphi(term, b + 1, 1, primes, &mut total)  ;
acc += muval1 as i64 * total ;
}
i+=2 ; }
acc
}
#[inline]
pub fn special_leaves_type_1( b : usize , interval : usize ,  m1 : &mut [usize] , n : usize, pp : usize , m : u64 ,
 	 interval_boundaries : &[usize], mu : &[isize], count : &mut i64, phi : &[u64], counter : &[i32]) { 
if m1[b] % 2 == 0 { m1[b] -= 1;} 
let criterion = n / pp ; 
while m1[b] > criterion {
   let y  = (m / (m1[b] as u64 * pp as u64 )) as usize  ; //print!("y = {} ",y) ;
   if y > interval_boundaries[interval + 1] - 2 { return ;} 
   let   muvalue = mu[((m1[b]+1) >> 1) as usize]; 
   if muvalue.abs() > pp as isize  {
        *count -=  
        muvalue.signum() as i64 * (phi[b] as i64 + cnt_query((y + 1 - interval_boundaries[interval]), counter) as i64);
   } 
   m1[b] -= 2 ;
      }
  } 
 	 
 #[inline]
pub fn special_leaves_type_2_initialize( b : usize, pb : usize, m : u64, t : &mut[usize], n : usize, pi : &[usize], a : usize, d2 : &mut[usize], count : &mut i64 )  { 
let    term = (m / (pb as u64 * pb as u64)) as usize;
t[b] = match  term {
   term if term <= pb =>  b + 2,
   term if term < n =>  pi[(term + 1) >> 1] + 1,
   _ => a + 1,
       };
d2[b] = t[b] - 1 ;
*count += (a - d2[b]) as i64 ;
}
#[inline]
pub fn hard(b :  usize , interval : usize, y : usize, interval_boundaries : &[usize], count : &mut i64, counter : &[i32], d2 : &mut [usize]) -> ( bool, u32) {
    if y + 1 >= interval_boundaries[interval+1] {return (true,0); }
    *count+= cnt_query((y + 1 - interval_boundaries[interval]),counter) as i64;
    d2[b] -= 1;
   ( false, 1 )
    }
   
 #[inline]  
pub fn easy_sparse(b :  usize , interval : usize, y : usize, n : usize, tt : &mut[u8], switch : &mut [bool],interval_boundaries : &[usize],
	 count : &mut i64, counter : &[i32], d2 : &mut [usize], pi : &[usize] ) -> bool  {
      if y >= n {
         if !switch[b] { switch[b]=true; return true; } else { tt[b] = 2 ; hard(b,interval,y,interval_boundaries,count,counter,d2); }
      }
      else {
         let l = pi[((y + 1) >> 1)] - b + 1;
         *count += l as i64;
         d2[b] -= 1;
      }
    false
    }
#[inline]   
pub fn easy_clustered(b :  usize , interval : usize, y : usize, n : usize, tt: &mut[u8], switch : &mut [bool], interval_boundaries : &[usize],
	 count : &mut i64, counter : &[i32], d2 : &mut [usize], m : u64, pi : &[usize], p : &[usize] ) -> bool  {
     if y >= n {
        if !switch[b] { switch[b]=true; return true; } else { tt[b] = 2 ; hard(b,interval,y,interval_boundaries,count,counter,d2); }
      }
     else
      {
     let   l = pi[((y + 1) >> 1)] - b + 1;
     let  term = m / (p[b + 1] as u64 * p[b + l] as u64);
     let  dprime = pi[((term + 1) >> 1) as usize];
     if p[dprime + 1] <= int_sqrt((m / p[b + 1] as u64) as usize) || dprime <= b  {
         tt[b] = 1;
         *count += l as i64;
         d2[b] -= 1 ;
        }
      else
    {
      *count +=  (l as u32 * (d2[b] - dprime) as u32) as i64 ;
      d2[b] = dprime;
      }
    }
    false
      } 
    
  #[inline]  
    pub fn special_leaves_type_2(b: usize, interval : usize, d2 : &mut[usize], m : u64, p : &[usize], tt : &mut[u8],
    	 n : usize, switch : &mut[bool], interval_boundaries : &[usize], count : &mut i64, counter : &[i32], pi : &[usize] )  -> u32 {
    let mut s2bprimes= 0;
     loop {
      if d2[b] == b + 1 {return s2bprimes; }
      else
      {
       let y = (m / (p[b + 1] as u64 * p[d2[b]] as u64)) as usize;
       match tt[b] {
          0 => { if easy_clustered(b, interval,y,n,tt,switch,interval_boundaries,count,counter,d2,m,pi,p) { return s2bprimes;} continue ; } ,
          1 => { if easy_sparse(b,interval,y,n,tt,switch,interval_boundaries,count,counter,d2,pi) { return s2bprimes; }  continue ; } ,
          _ => { let h = hard(b,interval,y,interval_boundaries,count,counter,d2) ;
                 if (interval > 0 || counter[1] > 0) && h.0 { return s2bprimes;} s2bprimes += h.1 ; continue ; } ,
}
      }             
      }
     } 

  #[inline]     
 pub fn sieve2( x : usize, y : usize,   p : &[usize],  block : &mut BitVec )
  { 
block.clear();
   let mut i : usize = 1; 
    while p[i] * p[i] <= y {
    let mut offset = modneg(1-x as i64,p[i]);
       while offset <= 1+(y-x)  { 
    	 block.set(offset,true);
         offset += p[i];
        }
       i+=1 ; 
    } 
}
 #[inline]
 pub fn p2(interval : usize, u : &mut usize, v :  &mut usize, n : usize, w : &mut usize,
   block : &mut BitVec , p : &[usize], m : u64 , interval_boundaries : &[usize], phi2 : &mut i64, counter : &[i32], a : usize  ) -> u32
      { 
  let mut p2primes = 0;
loop { 
    if *u <= n { *phi2 -= (*v as i64 * (*v - 1) as i64) >> 1 ; return p2primes;}
   	if *u  < *w  {
   		*w = cmp::max(2,*u -n); 
   		sieve2(*w,*u+1,p,block) ;
   		    } 
	let index = *u-*w+1;
    if !block[index] {
    let y  = (m / (*u as u64)) as usize;
    if y +1 >= interval_boundaries[interval + 1] { return p2primes; }  
    *phi2 += (cnt_query((y + 1 - interval_boundaries[interval]), counter) as usize + a) as i64 - 1;
    p2primes += 1; 
    *v += 1;
    }
    *u -= 2; 
}  
} 

