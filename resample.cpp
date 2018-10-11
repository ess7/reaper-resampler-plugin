/*
    Derived from WDL - resample.cpp
    Copyright (C) 2010 and later Cockos Incorporated
    This software is provided 'as-is', without any express or implied
    warranty.  In no event will the authors be held liable for any damages
    arising from the use of this software.
    Permission is granted to anyone to use this software for any purpose,
    including commercial applications, and to alter it and redistribute it
    freely, subject to the following restrictions:
    1. The origin of this software must not be misrepresented; you must not
       claim that you wrote the original software. If you use this software
       in a product, an acknowledgment in the product documentation would be
       appreciated but is not required.
    2. Altered source versions must be plainly marked as such, and must not be
       misrepresented as being the original software.
    3. This notice may not be removed or altered from any source distribution.
      
    You may also distribute this software under the LGPL v2 or later.
*/

#include <math.h>
#include <emmintrin.h>
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <stdio.h>

#include "resample.h"

#define PI 3.1415926535897932384626433832795

double i0(double);

void WDL_Resampler::BuildLowPass(const double filtpos, const double beta, const int interpsize) {
	const double i0betainv = 1.0/i0(beta);
	
	m_sincoversize = interpsize;
	const int wantsize=m_sincsize;
	const int wantinterp=m_sincoversize;
	if (m_filter_ratio!=filtpos || m_filter_coeffs_size != wantsize || m_lp_oversize != wantinterp) {
		m_lp_oversize = wantinterp;
		m_filter_ratio = filtpos;
		
		const int allocsize = wantsize*(m_lp_oversize + 3);
		WDL_SincFilterSample *cfout = m_filter_coeffs.Resize(allocsize);
		if (m_filter_coeffs.GetSize() == allocsize) {
			m_filter_coeffs_size = wantsize;
			
			const double dsincpos = PI*filtpos;
			const int hwantsize = wantsize/2 - 1;
			
			double filtpower = 0.0;
			WDL_SincFilterSample *ptrout = cfout;
			for (int slice = -1; slice <= wantinterp + 1; slice++) {
				const double frac = slice/(double)wantinterp;
				const int center_x = slice == 0 ? hwantsize+1 : slice == wantinterp ? hwantsize : -1;
				
				for (int x = 0; x < wantsize; x++) {
					if (x == center_x)  {
						*ptrout++ = 1.0;
					} else {
						const double t = x-1 + frac - hwantsize;
						const double sincpos = dsincpos * t;
						// kaiser*sinc
						const double tnorm = t/hwantsize;
						const double val = (tnorm >= -1.0 && tnorm <= 1.0)
							? (i0(beta*sqrt(1.0 - tnorm*tnorm))*i0betainv) * sin(sincpos)/sincpos
							: 0.0;
						if (slice >= 0 && slice < wantinterp) {
							filtpower+=val;
						}
						*ptrout++ = (WDL_SincFilterSample)val;
					}
				}
			}
			filtpower = wantinterp/(filtpower+1.0);
			for (int x = 0; x < allocsize; x ++) {
				cfout[x] = (WDL_SincFilterSample) (cfout[x]*filtpower);
			}
		} else {
			m_filter_coeffs_size = 0;
		}
	}
}

const WDL_SincFilterSample *WDL_Resampler::GetFilterCoeff(Cache *cache) {
	if (cache != NULL) {
		for (;;) {
			if (cache->state >= Cache::State::PROCESSING) {
				while (cache->state != Cache::State::READY) {
					Sleep(1);
				}
				m_filter_coeffs_size = cache->coeff != NULL ? m_sincsize : 0;
				m_lp_oversize = m_sincoversize = cache->interpsize;
				return cache->coeff;
			} else if (InterlockedCompareExchange((LONG *)&cache->state, Cache::State::PROCESSING, Cache::State::EMPTY)
					== Cache::State::EMPTY) {
				break;
			}
		}
	}
	
	#ifdef WDL_RESAMPLE_FULL_SINC_PRECISION
		const int atten = 180;
	#else
		const int atten = 160;
	#endif
	// from kaiserord
	const double beta = 0.1102*(atten-8.7);
	const double width = m_sratein*(atten-7.95)/(2.0*PI*2.285*m_sincsize);	
	#ifdef CUSTOM_FILTER
		// no aliasing cutoff
		const double fc = 0.5*((m_sratein <= m_srateout ? m_sratein : m_srateout) - width);
	#else
		// WDL cutoff
		const double fc = m_sratein <= m_srateout ? 0.5*m_sratein : (0.5/1.03)*m_srateout;
	#endif
	
	int interpsize;
	if (cache != NULL) {
		if (cache->i_out == SRATE_ARB) {
			interpsize = INTERPSIZE;
		} else if ((cache->i_in & 1) != (cache->i_out & 1)) {
			// n*44100 <=> m*48000
			interpsize = (int)(m_srateout*(1.0/300.0));  // m_srateout/gcd(44100, 48000)
		} else {
			// *n or /n
			interpsize = m_srateout > m_sratein ? (int)(m_srateout/m_sratein) : 1;
		}
		cache->interpsize = interpsize;
	} else {
		interpsize = INTERPSIZE;
	}
	BuildLowPass(fc/(0.5*m_sratein), beta, interpsize);
	char buf[128];
	sprintf(buf, "BuildLowPass %.2f->%.2f beta:%.2f fc:%.2f cache:(%d %d) size:%d interp:%d", m_sratein, m_srateout, beta, fc,
		cache != NULL ? cache->i_in : -1, cache != NULL ? cache->i_out : -1, m_filter_coeffs_size, m_lp_oversize);
	OutputDebugString(buf);
	
	const WDL_SincFilterSample *cfout = m_filter_coeffs.Get();
	if (cache != NULL) {
		const int size = sizeof(WDL_SincFilterSample)*m_filter_coeffs_size*(m_lp_oversize + 3);
		cache->coeff = (WDL_SincFilterSample *)malloc(size);
		if (cache->coeff != NULL) {
			memcpy(cache->coeff, cfout, size);
		} else {
			m_filter_coeffs_size = 0;
		}
		cache->state = Cache::State::READY;
		return cache->coeff;
	} else {
		return cfout;
	}
}

inline void WDL_Resampler::SincSampleZOH2N(WDL_ResampleSample *outptr, const WDL_ResampleSample *inptr, double fracpos, int nch, const WDL_SincFilterSample *filter, int filtsz) {
	const int oversize=m_lp_oversize;
	fracpos *= oversize;
	const int ifpos = lround(fracpos);
	const WDL_SincFilterSample *fptr = filter + (oversize-ifpos+1) * filtsz;
	
	const int blksz = nch == 2 ? 4 : nch == 4 ? 2 : 1;
	
	__m128d sum[WDL_RESAMPLE_MAX_NCH/2];
	for (int ch = 0; ch < nch/2; ch++) {
		sum[ch] = _mm_setzero_pd();
	}
	// assumes filtsz % blksz == 0
	while (filtsz) {
		for (int i = 0; i < blksz; i++) {
			const double _f = *fptr;
			const __m128d f = _mm_load1_pd(&_f);
			for (int ch = 0; ch < nch/2; ch++) {
				sum[ch] = _mm_add_pd(sum[ch], _mm_mul_pd(_mm_loadu_pd(inptr + 2*ch), f));
			}
			inptr += nch;
			fptr++;
		}
		filtsz -= blksz;
	}
	for (int ch = 0; ch < nch/2; ch++) {
		_mm_storeu_pd(outptr + 2*ch, sum[ch]);
	}
}

template<typename T, typename V>
static inline V quadinterp(T x, V y1, V y2, V y3) {
	// 0.5*y1*x*(x-1.0) - y2*(x+1.0)*(x-1.0) + 0.5*y3*(x+1.0)*(x);
	const T xm1 = x - 1.0;
	const T xp1 = x + 1.0;
	return 0.5*x*(y1*xm1 + y3*xp1) - y2*xp1*xm1;
}

inline void WDL_Resampler::SincSampleQuad2N(WDL_ResampleSample *outptr, const WDL_ResampleSample *inptr, double fracpos, int nch, const WDL_SincFilterSample *filter, int filtsz) {
	const int oversize=m_lp_oversize;
	fracpos *= oversize;
	const int ifpos = lround(fracpos);
	fracpos -= ifpos;
	for (int ch = 0; ch < nch; ch += 2) {
		const WDL_SincFilterSample *fptr1 = filter + (oversize-ifpos+2) * filtsz;
		const WDL_SincFilterSample *fptr2 = fptr1 - filtsz;
		const WDL_SincFilterSample *fptr3 = fptr2 - filtsz;
		const WDL_ResampleSample *iptr = inptr + ch;
		__m128d sum1 = _mm_setzero_pd();
		__m128d sum2 = _mm_setzero_pd();
		__m128d sum3 = _mm_setzero_pd();
		int i = filtsz;
		while (i--) {
			const __m128d inp = _mm_loadu_pd(iptr);
			const double f[3] = {*fptr1, *fptr2, *fptr3};
			sum1 = _mm_add_pd(sum1, _mm_mul_pd(inp, _mm_load1_pd(f  )));
			sum2 = _mm_add_pd(sum2, _mm_mul_pd(inp, _mm_load1_pd(f+1)));
			sum3 = _mm_add_pd(sum3, _mm_mul_pd(inp, _mm_load1_pd(f+2)));
			iptr += nch;
			fptr1++;
			fptr2++;
			fptr3++;
		}
		double out[6];
		_mm_storeu_pd(out  , sum1);
		_mm_storeu_pd(out+2, sum2);
		_mm_storeu_pd(out+4, sum3);
		outptr[ch]   = quadinterp(fracpos, out[0], out[2], out[4]);
		outptr[ch+1] = quadinterp(fracpos, out[1], out[3], out[5]);
	}
}

int WDL_Resampler::ResampleOut(WDL_ResampleSample *out, int nsamples_in, int nsamples_out, int nch) 
{
  if (nch > WDL_RESAMPLE_MAX_NCH || nch < 1) return 0;
  if (m_filtercnt>0) return 0;  // not implemented
  
  // prevent the caller from corrupting the internal state
  m_samples_in_rsinbuf += nsamples_in < m_last_requested ? nsamples_in : m_last_requested; 

  int rsinbuf_availtemp = m_samples_in_rsinbuf;

  if (nsamples_in < m_last_requested) // flush out to ensure we can deliver
  {
    int fsize=(m_last_requested-nsamples_in)*2 + m_sincsize*2;

    int alloc_size=(m_samples_in_rsinbuf+fsize)*nch;
    WDL_ResampleSample *zb=m_rsinbuf.Resize(alloc_size,false);
    if (m_rsinbuf.GetSize() == alloc_size)
    {
      memset(zb+m_samples_in_rsinbuf*nch,0,fsize*nch*sizeof(WDL_ResampleSample));
      rsinbuf_availtemp = m_samples_in_rsinbuf+fsize;
    }
  }

  int ret=0;
  double srcpos=m_fracpos;
  double drspos = m_ratio;
  WDL_ResampleSample *localin = m_rsinbuf.Get();

  WDL_ResampleSample *outptr=out;

  int ns=nsamples_out;

  int outlatadj=0;

  if (m_sincsize) // sinc interpolating
  {
	Cache *cache = GetCache();
    const WDL_SincFilterSample *filter = GetFilterCoeff(cache);
    int filtsz = m_filter_coeffs_size;
    int filtlen = rsinbuf_availtemp - filtsz;
    outlatadj = filtsz/2-1;
	const bool useZOH = cache != NULL ? cache->i_out != SRATE_ARB : false;

	#define SINCSAMPLE_LOOP(SINCSAMPLE, NCH) \
	  while (ns--) { \
        int ipos = (int)srcpos; \
        if (ipos >= filtlen-1) { break; } \
		SINCSAMPLE(outptr, localin + ipos*NCH, srcpos-ipos, NCH, filter, filtsz); \
        outptr += NCH; \
        srcpos += drspos; \
        ret++; \
      }

    if (nch % 2 == 0)
    {
      if (useZOH) {
	    switch (nch) {  // common channel counts
		  case 2:
		    SINCSAMPLE_LOOP(SincSampleZOH2N, 2)
		    break;
		  case 6:
		    SINCSAMPLE_LOOP(SincSampleZOH2N, 6)
		    break;
		  case 8:
		    SINCSAMPLE_LOOP(SincSampleZOH2N, 8)
		    break;
		  default:
		    SINCSAMPLE_LOOP(SincSampleZOH2N, nch)
		}
	  } else {
	    SINCSAMPLE_LOOP(SincSampleQuad2N, nch)
	  }
    }
	else {
		// not yet implemented
		return 0;
	}
	if (useZOH) {
		// snap srcpos to 1/m_lp_oversize grid
		srcpos = round(srcpos*m_lp_oversize)/m_lp_oversize;
	}
  }
  else
  {
	// not implemented
    return 0;
  }
  
  if (ret>0 && rsinbuf_availtemp>m_samples_in_rsinbuf) // we had to pad!!
  {
    // check for the case where rsinbuf_availtemp>m_samples_in_rsinbuf, decrease ret down to actual valid samples
    double adj=(srcpos-m_samples_in_rsinbuf + outlatadj) / drspos;
    if (adj>0)
    {
      ret -= (int) (adj + 0.5);
      if (ret<0)ret=0;
    }
  }

  int isrcpos=(int)srcpos;
  if (isrcpos > m_samples_in_rsinbuf) isrcpos=m_samples_in_rsinbuf;
  m_fracpos = srcpos - isrcpos;
  m_samples_in_rsinbuf -= isrcpos;
  if (m_samples_in_rsinbuf <= 0) m_samples_in_rsinbuf=0;
  else
    memmove(localin, localin + isrcpos*nch,m_samples_in_rsinbuf*sizeof(WDL_ResampleSample)*nch);

  return ret;
}

inline WDL_Resampler::SampleRate WDL_Resampler::GetSrateIndex(float srate) {
	if (srate == 44100.0f) {
		return SRATE_44p1;
	} else if (srate == 48000.0f) {
		return SRATE_48;
	} else if (srate == 2*48000.0f) {
		return SRATE_96;
	} else if (srate == 2*44100.0f) {
		return SRATE_88p2;
	} else if (srate == 4*44100.0f) {
		return SRATE_176p4;
	} else if (srate == 4*48000.0f) {
		return SRATE_192;
	} else {
		return SRATE_UNK;
	}
}

WDL_Resampler::Cache *WDL_Resampler::GetCache() {
	SampleRate i_in, i_out;
	i_in = GetSrateIndex(m_sratein);
	if (i_in == SRATE_UNK) {
		return NULL;
	}
	i_out = GetSrateIndex(m_srateout);
	if (i_out == SRATE_UNK) {
		if (m_srateout >= m_sratein) {
			i_out = SRATE_ARB;
		} else {
			return NULL;
		}
	}
	
	Cache *cache;
	switch (m_sincsize) {
		case 192:
			cache = &m_cache[SIZE_192][i_in][i_out];
			break;
		case 384:
			cache = &m_cache[SIZE_384][i_in][i_out];
			break;
		case 64:
			cache = &m_cache[SIZE_64 ][i_in][i_out];
			break;
		case 512:
			cache = &m_cache[SIZE_512][i_in][i_out];
			break;
		default:
			return NULL;
	}
	cache->i_in  = i_in;
	cache->i_out = i_out;
	return cache;
}

decltype(WDL_HeapBuf::pResize) WDL_HeapBuf::pResize;
decltype(WDL_Resampler::pResampleOut) WDL_Resampler::pResampleOut(&WDL_Resampler::ResampleOut);
decltype(WDL_Resampler::m_cache) WDL_Resampler::m_cache;
