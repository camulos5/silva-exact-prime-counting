#![feature(negate_unsigned,box_syntax,step_by)]
#![feature(alloc_jemalloc)]

extern crate alloc_jemalloc;

const SIGNBIT : i32 = 1<<31;

pub fn modneg( x : i64, y : usize) -> usize {
	(( x % y as i64) + y as i64) as usize % y
	}

 pub fn signum(n : isize) -> i8 {
 	if n < 0 { -1 }
 	else if n > 0 { 1 }
 	else { 0 }
 	}
 
 
 pub fn nabs(n : isize) -> isize {
 	if n < 0 { -n } else { n }
 	//(n & !(1<<31)) as u32
 }
 
 pub fn maxm(m : isize , n :isize) -> isize {
 	if m>= n {m} else {n}
 }


pub fn int_sqrt(n :usize) -> usize {
	((n as f64).sqrt()).floor() as usize
}  


pub fn cnt_query(mut pos : usize, counter : &[i32]) -> u32 {
let mut  acc : i64 = counter[pos] as i64 ; //& !SIGNBIT) as i64; // as isize;
pos+=1 ;
pos &= pos -1 ;
while pos > 1 {
acc += counter[pos - 1] as i64;
//println!("acc = {}",acc);//& !SIGNBIT ;
pos &= pos - 1 ; }
(acc as i32 & !SIGNBIT) as u32
}


pub fn cnt_update(mut pos : usize, counter : &mut[i32], interval_length : usize ) {
    counter[pos] ^= SIGNBIT;
    while pos < interval_length  {
      counter[pos]-=1;
      pos |= pos + 1; 
    }
    }


pub fn cnt_init(counter: &mut[i32], interval_length : usize) {
for i in 0..interval_length {
      counter[i] = (i as i32 +1)  & !(i as i32);
      }
}


pub fn interval_clear( b : usize , offset : &mut[usize], counter : &mut[i32], interval_length : usize, pp : usize, lc :  u8)  {
   let mut j  = offset[b] ;
// for x in  (offset[b] as usize..interval_length).step_by(pp as usize)  {
// 	if counter[x] >0 {cnt_update(x,counter,interval_length) ; } 
// 	} 
   while j < interval_length {
   	if counter[j] > 0 { cnt_update( j, counter, interval_length ) ; }
   	j += pp ;
   	}
   offset[b] = modneg( offset[b] as i64 - (1 << lc) , pp );
   }


 pub fn rphi(n: i64,  mut a : usize, bit : i8 , pp: &[usize],  total : &mut i64 )  {
	loop {
		if a==1 {
			*total +=  bit as i64 * n ; 
			  break;
			}
			else if n < pp[a] as i64{
			*total += bit as i64 ;
				break;
				}
			else {
				a-=1;
				rphi(n / (pp[a] as i64), a, -bit, pp, total) ;
					}
				}
		}
 
 /*
 
 pub fn foo(x: Box<Vec<bool>>) -> Vec<bool> {
    *x
}
 
 
 pub fn bar( y:  Box<Vec<i8>>) -> Vec<i8> {
    *y
}
 
  
 pub fn baz(z: Box<&[u32]>) -> &[u32] {
    *z
}
*/
 
 
 pub fn special_leaves_type_1( b : usize , interval : usize ,  m1 : &mut [usize] , n : usize, pp : usize , m : u64 ,
 	 interval_boundaries : &[usize], mu : &[isize], count : &mut i64, phi : &[u64], counter : &[i32]) { 
    if m1[b] % 2 ==0 { m1[b]-=1;} //m1[b] += m1[b] % 2 - 1;
  let criterion = n / pp ; //print!("m1[b] = {}, criterion = {}, interval_boundaries[interval] = {} ",m1[b], criterion, interval_boundaries[interval]) ;
      while m1[b] > criterion {
   let    y  = (m / (m1[b] as u64 * pp as u64 )) as usize  ; //print!("y = {} ",y) ;
    if y > interval_boundaries[interval + 1] - 2 { return ;} 
 let   muvalue = mu[((m1[b]+1) >> 1) as usize]; 
    if nabs(muvalue) > pp as isize{
        *count -=  
    signum(muvalue) as i64 * (phi[b] as i64 + cnt_query((y + 1 - interval_boundaries[interval]), counter) as i64);
   } 
     m1[b] -= 2 ;
      }
  } 
 	 
 
 pub fn sieve2( x : usize, y : usize,   p : &[usize],  block : &mut [bool] )
  { 
for b in 0..block.len() { block[b] = false ; } 
// for _ in block.to_vec() { false;}
//block = &mut [false;block.len()] ;
// block =  block.iter_mut().map(move |_| false).collect::<&[bool]>();// { b = &mut false;}
//println!("{:?}",block);
   let mut i : usize = 1; 
    while p[i] * p[i] <= y {
    let mut offset = modneg(1-x as i64,p[i]);
//    (offset..2+(y-x) as usize).step_by(p[i] as usize).map(|z| blockc[z] = true).collect::<Vec<_>>() ;
      while offset <= 1+(y-x)  { 
	 block[offset] = true;
        offset += p[i];
      }
   i+=1 ; 
    } 
}
 
 pub fn p2(interval : usize, p2primes :  &mut usize, u : &mut usize, v :  &mut usize, n : usize, w : &mut usize,
   block : &mut [bool] , p : &[usize], m : u64 , interval_boundaries : &[usize], phi2 : &mut i64, counter : &[i32], a : usize  )
      { 
  *p2primes = 0;
loop { 
    if *u <= n { *phi2 -= (*v as i64 * (*v - 1) as i64) >> 1 ; return;}
   	if *u  < *w  {
   		*w = maxm(2,(*u -n) as isize) as usize; 
   		    sieve2(*w,*u+1,p,block) ;
   		    } 
	let index = *u-*w+1;
    if !block[index] {// println!("*u = {}", *u) ;
    let y  = (m / (*u as u64)) as usize;
    if y +1 >= interval_boundaries[interval + 1] { return ; }  
    *phi2+= (cnt_query((y + 1 - interval_boundaries[interval]), counter) as usize + a) as i64 - 1;
    *p2primes+=1; 
    *v+=1;
    }
    *u-=2; 
}  
} 

 
 pub fn initialize_arrays(ll : usize, mu : &mut[isize], pix : &mut usize, pi : &mut[usize], p : &mut[usize])  {
//mu[0] = 0;
//for i in 1..mu.len() { mu[i] =1 ; }
 	for j in  2..mu.len() {
 		if mu[j] == 1 {
 			for i in (j..ll+2).step_by(2 * j - 1) {
 			 mu[i] = match mu[i] {
 				1 => 1 - 2 * j as isize,
 				_ => -mu[i] ,
 			} ;	
 			}
 		}
 	}
 	for j in 2..int_sqrt(ll<<1) {
 		if mu[j] == 1- 2*j as isize{
 			for i in ((2*j*j -2*j +1)..ll+2).step_by( 4*j*j -4*j +1 ) {
 			mu[i] = 0 ;	
 			}
 		}
 	}
 pi[1]=0;pi[2]=2;
 p[0]=1;p[1]=2;
 *pix = 1 ;
let mut i = 2 ; for j in mu.iter().skip(2) {
// for i in 2..mu.len() {
if *j == 1 - 2* i {
// if mu[i] == 1 - 2*i as isize {
 	*pix+=1;
 	p[*pix] = 2*i as usize -1 ;
}
 pi[i as usize]=*pix;	
i+=1 ;
 }
 } 
 

 pub fn ordinary_leaves(n : usize, count : &mut i64, mu : &[isize], m : u64) {
  for i in (1..(n + 1)).step_by(2) {
  	*count+= signum(mu[(i+1) >> 1]) as i64 * (m as i64/ i as i64) ;
  }
for i in (2..(n + 1)).step_by(4) {
		*count -= (m as i64 / i as i64) * signum(mu[((i >> 1)  + 1) >> 1]) as i64 ;
	}
}       


pub fn special_leaves_type_1_substitute(b : usize, p : &[usize], n : usize, mu : &[isize], m : u64, count : &mut i64) {
let pp = p[b + 1]    ;
let mut j = maxm((n / pp)as isize, pp as isize) as usize +1 ;
if j % 2 == 0 { j+=1;} 
for i in (j..(n+1)).step_by(2) {
let muval1 = signum(mu[ ( i+1) >>1 ] ) ;
if nabs(mu[(i + 1) >> 1]) > pp as isize && muval1 != 0 {
let term = (m / (i * pp) as u64) as i64;
let mut total : i64 = 0;
rphi(term, b + 1, 1, p, &mut total)  ;
*count-= muval1 as i64 * total;
}
}
}


pub fn special_leaves_type_2_initialize( b : usize, pb : usize, m : u64, t : &mut[usize], n : usize, pi : &[usize], a : usize, d2 : &mut[usize], count : &mut i64 )  { 
let    term = (m / (pb as u64 * pb as u64)) as usize;
 t[b] = match  term {
   term if term <= pb =>  b + 2,
    term if term < n =>  pi[(term + 1) >> 1] + 1,
    _ => a + 1,
       };
    d2[b] = t[b] - 1 ;
    *count+= (a - d2[b]) as i64 ;
}

  pub fn hard(b :  usize , interval : usize, y : usize, interval_boundaries : &[usize], count : &mut i64, counter : &[i32], s2bprimes : &mut usize, d2 : &mut [usize]) -> bool {
      if y + 1 >= interval_boundaries[interval+1] {return true; }
       *count+= cnt_query((y + 1 - interval_boundaries[interval]),counter) as i64;
       *s2bprimes += 1 ;
    d2[b] -= 1;
    false
    }
   
   
   pub fn easy_sparse(b :  usize , interval : usize, y : usize, n : usize, tt : &mut[u8], switch : &mut [bool],interval_boundaries : &[usize], count : &mut i64, counter : &[i32], s2bprimes : &mut usize, d2 : &mut [usize], pi : &[usize] ) -> bool  {
      if y >= n {
        if !switch[b] { switch[b]=true; return true; } else { tt[b] = 2 ; hard(b,interval,y,interval_boundaries,count,counter,s2bprimes,d2); }
      }
      else {
     let   l = pi[((y + 1) >> 1)] - b + 1;
     *count += l as i64;
    d2[b] -= 1;
    }
    false
    }   
  
   
   pub fn easy_clustered(b :  usize , interval : usize, y : usize, n : usize, tt: &mut[u8], switch : &mut [bool], interval_boundaries : &[usize], count : &mut i64, counter : &[i32], s2bprimes : &mut usize, d2 : &mut [usize], m : u64, pi : &[usize], p : &[usize] ) -> bool  {
     if y >= n {
        if !switch[b] { switch[b]=true; return true; } else { tt[b] = 2 ; hard(b,interval,y,interval_boundaries,count,counter,s2bprimes,d2); }
      }
      else
      {
     let   l = pi[((y + 1) >> 1)] - b + 1;
      let  term2 = m / (p[b + 1] as u64 * p[b + l] as u64);
      let  dprime = pi[((term2 + 1) >> 1) as usize];
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
    
    
    pub fn special_leaves_type_2(b: usize, interval : usize, s2bprimes : &mut usize, d2 : &mut[usize], m : u64, p : &[usize], tt : &mut[u8], n : usize, switch : &mut[bool], interval_boundaries : &[usize], count : &mut i64, counter : &[i32], pi : &[usize] )   {
    *s2bprimes= 0;
     loop {
      if d2[b] == b + 1 {return ; }
      else
      {
       let y = (m / (p[b + 1] as u64 * p[d2[b]] as u64)) as usize;
       match tt[b] {
          0 => { if easy_clustered(b, interval,y,n,tt,switch,interval_boundaries,count,counter,s2bprimes,d2,m,pi,p) { return;} continue ; } ,
          1 => { if easy_sparse(b,interval,y,n,tt,switch,interval_boundaries,count,counter,s2bprimes,d2,pi) { return ; }  continue ; } ,
          _ => { if hard(b,interval,y,interval_boundaries,count,counter,s2bprimes,d2) {return;} continue ; } ,
          }
      }                           
   }
     } 
      