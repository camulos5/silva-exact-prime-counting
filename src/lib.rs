#![feature(negate_unsigned,box_syntax,step_by)]
#![feature(alloc_jemalloc)]

extern crate alloc_jemalloc;

const SIGNBIT : i32 = 1<<31;

pub fn cnt_query(mut pos : usize, counter : &[i32]) -> u32 {
let mut  result=counter[pos] & !SIGNBIT ; // as isize;
pos+=1 ;
pos &= pos -1 ;
while pos > 1 {
result+=counter[pos - 1] & !SIGNBIT ;
pos &= pos - 1 ; }
(result & !SIGNBIT) as u32
}


pub fn cnt_update(mut pos : usize, counter : &mut[i32], interval_length : usize ) {
    counter[pos] |= 1<<31;
    while pos < interval_length  {
      counter[pos]-=1;
      pos |= pos + 1; 
    }
    }


pub fn cnt_init(counter: &mut[i32], interval_length : usize) {
for i in 0..interval_length {
      counter[i] = ((i as i32 +1)  & !(i as i32)) as i32 ;
      }
}


pub fn modneg( x : i64, y : usize) -> usize {
	(( x % y as i64) + y as i64) as usize % y
	}


pub extern "system" fn intervalclear( b : usize , offset : &mut[u32], counter : &mut[i32], interval_length : usize, pp : u32, lc :  u8)  {
   let mut j  = offset[b] as usize;
// for x in  (offset[b] as usize..interval_length).step_by(pp as usize)  {
// 	if counter[x] >0 {cnt_update(x,counter,interval_length) ; } 
// 	} 
   while j < interval_length {
   	if counter[j] > 0 { cnt_update( j, counter, interval_length ) ; }
   	j += pp as usize ;
   	}
   offset[b] = modneg( offset[b] as i64 - (1 << lc) , pp as usize) as u32;
   }


 pub fn rphi(n: i64,  mut a : usize, bit : i8 , pp: &[u32],  total : &mut i64 )  {
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
 
 
 pub fn signum(n : i32) -> i8 {
 	if n < 0 { -1 }
 	else if n > 0 { 1 }
 	else { 0 }
 	}
 
 
 pub fn nabs(n : i32) -> i32 {
 	if n < 0 { -n } else { n }
 	//(n & !(1<<31)) as u32
 }
 
 pub fn maxm(m : i32 , n :i32) -> i32 {
 	if m>= n {m} else {n}
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
 
 
 pub fn s1b( b : usize , interval : usize ,  m1 : &mut [u32] , n : u32, pp : u32 , m : u64 ,
 	 interval_boundaries : &[u32], mu : &[i32], count : &mut i64, phi : &[u64], counter : &[i32]) { 
    if m1[b] % 2 ==0 { m1[b]-=1;} //m1[b] += m1[b] % 2 - 1;
  let criterion : u32 = n / pp ; //print!("m1[b] = {}, criterion = {}, interval_boundaries[interval] = {} ",m1[b], criterion, interval_boundaries[interval]) ;
      while m1[b] > criterion {
   let    y : u32 = (m / (m1[b] as u64 * pp as u64 )) as u32  ; //print!("y = {} ",y) ;
    if y > interval_boundaries[interval + 1] - 2 { return ;} 
 let   muvalue : i32 = mu[((m1[b]+1) >> 1) as usize]; 
    if nabs(muvalue) as u32 > pp {
        *count -=  
    signum(muvalue) as i64 * (phi[b] as i64 + cnt_query((y + 1 - interval_boundaries[interval]) as usize,counter) as i64);
   } 
     m1[b] -= 2 ;
      }
  } 
 	 
 
 pub fn sieve2( x : u32, y : u32,   p : &[u32],   blockc : &mut [bool] )
  { //println!("sieve2");
// let  intervl : usize = (y - x ) as usize;
// println!("blockc in sieve 2 = {:?}, blockc.len() = {}",blockc, blockc.len()) ;
//let mut block = blockc ;
for b in 0..blockc.len() { blockc[b] = false ; }
//for _ in blockc.iter() { false; } // compiles but doesn't initialize to false
//  println!("blockc in sieve 2 after initialization = {:?}, blockc.len() = {}",blockc, blockc.len()) ;
   let mut i : usize = 1; 
    while p[i] * p[i] <= y  {
    let mut offset = modneg(1-x as i64,p[i] as usize);
//    (offset..2+(y-x) as usize).step_by(p[i] as usize).map(|z| blockc[z] = true).collect::<Vec<_>>() ;
      while offset <= 1+(y-x) as usize  { 
	 blockc[offset] = true;
        offset += p[i] as usize;
      }
   i+=1 ; 
    } 
}
 // block doesn't need to be global 
 
 pub fn p2(interval : usize, p2primes :  &mut u32, u : &mut u32, v :  &mut usize, n : u32, w : &mut u32,
   block : &mut [bool] , p : &[u32], m : u64 , interval_boundaries : &[u32], phi2 : &mut i64, counter : &[i32], a : usize  )
      { 
  *p2primes = 0;
loop { 
    if *u <= n { *phi2 -= (*v as i64 * (*v - 1) as i64) >> 1 ; return;}
   	if *u  < *w  {
   		*w = maxm(2,(*u-n) as i32)as u32; 
//   		for mut b in block.iter() { b = &false ;}
   		    sieve2(*w,*u+1,p,block) ;
   		    } //*println!("block upon return is {:?}",block);*/}
	let index = (*u-*w+1) as usize ;
    if !block[index] {// println!("*u = {}", *u) ;
    let y : u32 = (m / (*u as u64)) as u32;
    if y +1 >= interval_boundaries[interval + 1] { return ; }  
    *phi2+= (cnt_query((y + 1 - interval_boundaries[interval]) as usize,counter) as usize + a) as i64 - 1;
    *p2primes+=1; 
    *v+=1;
    }
    *u-=2; 
}  
} 


pub fn int_sqrt(n :usize) -> usize {
	((n as f64).sqrt()).floor() as usize
}  

 pub fn exp(n: u64, mut e:usize) -> u64 {
 	let mut num = n ;
 	while e > 1 {
 		num*=n;
 		e-=1;
 	}
 num
 }
 
 pub fn init_arrays(ll : usize, mu : &mut[i32], pix : &mut u32, pi : &mut[u32], p : &mut[u32])  {
for i in 1..mu.len() { mu[i] =1 ; }
 	for j in  2..mu.len() {
 		if mu[j] == 1 {
 			for i in (j..ll+2).step_by(2 * j - 1) {
 			 mu[i] = match mu[i] {
 				1 => 1 - 2 * j as i32,
 				_ => -mu[i] ,
 			} ;	
 			}
 		}
 	}
 	for j in 2..int_sqrt(ll<<1) {
 		if mu[j] == 1- 2*j as i32{
 			for i in ((2*j*j -2*j +1)..ll+2).step_by( 4*j*j -4*j +1 ) {
 			mu[i] = 0 ;	
 			}
 		}
 	}
 pi[1]=0;pi[2]=2;
 p[0]=1;p[1]=2;
 *pix = 1 ;
 for i in 2..mu.len() {
 if mu[i] == 1 - 2*i as i32 {
 	*pix+=1;
 	p[*pix as usize]=2*i as u32-1 ;
 }
 pi[i]=*pix;	
 }
 } 
 

 pub fn ord_leaves(n : u32, count : &mut i64, mu : &[i32], m : u64) {
  for i in (1..(n as usize +1)).step_by(2) {
  	*count+= signum(mu[(i+1) >> 1]) as i64 * (m as i64/ i as i64) ;
  }
for i in (2..(n as usize +1)).step_by(4) {
		*count -= (m as i64 / i as i64) * signum(mu[((i >> 1)  + 1) >> 1]) as i64 ;
	}
}       


pub fn s1b_subst(b : usize, p : &[u32], n : u32, mu : &[i32], m : u64, count : &mut i64) {
//PROCEDURE s1b_subst (b : cardinal) ;
//more general algorithm
//var term : int64 ; i, pp : cardinal ;   muval1 : shortint ;  total : int64 ;
//begin
let pp = p[b + 1]    ;
let mut j = maxm((n / pp)as i32, pp as i32) as usize +1 ;
if j % 2 == 0 { j+=1;} //j -= (j % 2) - 1 ;
for i in (j..(n+1)as usize).step_by(2) {
//while i <= n do  begin
let muval1 = signum(mu[ ( i+1) >>1 ] ) ;
if nabs(mu[(i + 1) >> 1]) as u32 > pp && muval1 != 0 {
let term = (m / (i as u32 * pp) as u64) as i64;
let mut total : i64 = 0;
rphi(term, b + 1, 1, p, &mut total)  ;
*count-= muval1 as i64 * total;
}
}
}


pub fn s2b_init( b : usize, pb : u32, m : u64, t : &mut[u32], n : u32, pi : &[u32], a : usize, d2 : &mut[u32], count : &mut i64, tt : &mut[u8] )  { 
//PROCEDURE s2b_init( b: integer);
//  var
//    term, pb: cardinal;
//  begin
//    pb := p[succ(b)];
let    term = (m / (pb as u64 * pb as u64)) as u32;
    if term <= pb {
      t[b] = b as u32 + 2 ; }
    else if term < n {
      t[b] = pi[(term as usize + 1) >> 1] + 1 ; }
    else   {
      t[b] = a as u32 + 1 ; }
    d2[b] = t[b] - 1 ;
    *count+= (a as u32 - d2[b]) as i64 ;
    tt[b] = 0;
}

  pub fn hard(b :  usize , interval : usize, y : u64, interval_boundaries : &[u32], count : &mut i64, counter : &[i32], s2bprimes : &mut u32, d2 : &mut [u32]) -> bool {
      if y + 1 >= interval_boundaries[interval+1] as u64 {return true; }
       *count+= cnt_query((y + 1 - interval_boundaries[interval]as u64) as usize,counter) as i64;
       *s2bprimes+=1 ;
    d2[b]-=1;
    false
    }
   
   
   pub fn easy_sparse(b :  usize , interval : usize, y : u64, n : u32, tt : &mut[u8], switch : &mut [bool],interval_boundaries : &[u32], count : &mut i64, counter : &[i32], s2bprimes : &mut u32, d2 : &mut [u32], pi : &[u32] ) -> bool  {
//  var l : cardinal ;
//  begin
//    easy_sparse:=false ;
      if y >= n as u64 {
        tt[b] = 2;
        if !switch[b] { switch[b]=true; return true; } else {hard(b,interval,y,interval_boundaries,count,counter,s2bprimes,d2); }
      }
      else {
     let   l = pi[((y + 1) >> 1) as usize] - b as u32+ 1;
     *count += l as i64;
    d2[b]-=1;
    }
    false
    }   
  
   
   pub fn easy_clustered(b :  usize , interval : usize, y : u64, n : u32, tt: &mut[u8], switch : &mut [bool], interval_boundaries : &[u32], count : &mut i64, counter : &[i32], s2bprimes : &mut u32, d2 : &mut [u32], m : u64, pi : &[u32], p : &[u32] ) -> bool  {
     if y >= n as u64{
//        tt[b] = 2;
        if !switch[b] { switch[b]=true; return true; } else { tt[b] = 2 ; hard(b,interval,y,interval_boundaries,count,counter,s2bprimes,d2); }
      }
      else
      {
     let   l = pi[((y + 1) >> 1) as usize] as usize - b + 1;
      let  term2 = m / (p[b + 1] as u64 * p[b + l] as u64);
      let  dprime = pi[((term2 + 1) >> 1) as usize];
          if p[dprime as usize + 1] <= int_sqrt((m / p[b + 1] as u64) as usize) as u32 || dprime <= b as u32 {
      tt[b] = 1;
      *count += l as i64;
    d2[b]-=1 ;
}
              else
    {
      *count +=  (l as u32 * (d2[b] - dprime)) as i64 ;
      d2[b] = dprime;
      }
    }
    false
      } 
    
    
    pub fn s2b_structured(b: usize, interval : usize, s2bprimes : &mut u32, d2 : &mut[u32], m : u64, p : &[u32], tt : &mut[u8], n : u32, switch : &mut[bool], interval_boundaries : &[u32], count : &mut i64, counter : &[i32], pi : &[u32] )   {
    *s2bprimes= 0;
     loop {
      if d2[b] == b as u32 + 1 {return ; }
      else
      {
       let y = m / (p[b + 1] as u64 * p[d2[b]as usize] as u64);
       match tt[b] {
         //0: begin if easy_clustered(b,interval,y)then exit; continue; end ;
          0 => { if easy_clustered(b, interval,y,n,tt,switch,interval_boundaries,count,counter,s2bprimes,d2,m,pi,p) { return;} continue ; } ,
          1 => { if easy_sparse(b,interval,y,n,tt,switch,interval_boundaries,count,counter,s2bprimes,d2,pi) { return ; }  continue ; } ,
          2 => { if hard(b,interval,y,interval_boundaries,count,counter,s2bprimes,d2) {return;} continue ; } ,
          _ => println!("match error"),
          }
      }                           
   }
     }          