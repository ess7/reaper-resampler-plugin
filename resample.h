/*
    Derived from WDL - resample.h
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

#ifndef _WDL_RESAMPLE_H_
#define _WDL_RESAMPLE_H_

#include <type_traits>
#include <stddef.h>

#include "config.h"

#ifdef WDL_RESAMPLE_FULL_SINC_PRECISION
	#define INTERPSIZE 1280
	typedef double WDL_SincFilterSample;
#else
	#define INTERPSIZE 512
	typedef float WDL_SincFilterSample;
#endif

typedef double WDL_ResampleSample;

#define WDL_RESAMPLE_MAX_NCH 64


template<typename>
class MemberFuncPtr;

template<typename Ret, typename Cls, typename... Args>
class MemberFuncPtr<Ret (Cls::*)(Args...)> {
using T = Ret (Cls::*)(Args...);
static_assert(!std::is_polymorphic<Cls>::value, "MemberFuncPtr does not support polymorphic types");
public:
	MemberFuncPtr() : funcPtr(NULL) {}
    MemberFuncPtr(T funcPtr) : funcPtr(funcPtr) {}
	MemberFuncPtr(void *ptr) {
        raw.ptr = ptr;
		raw.adj = 0;
	}
    void *getAddr() const {
        return this->raw.ptr;
    }
	Ret invoke(Cls *that, Args... args) const {
        return (that->*(funcPtr))(args...);
    }
private:
    union {
        struct {
            void *ptr;
            ptrdiff_t adj;
        } raw;
        T funcPtr;
    };
};


class WDL_HeapBuf {
public:
	void *Get() const { return m_size ? m_buf : NULL; } // returns NULL if size is 0
	int GetSize() const { return m_size; }
	void *Resize(int newsize, bool resizedown) {
		return pResize.invoke(this, newsize, resizedown);
	}
	static MemberFuncPtr<decltype(&WDL_HeapBuf::Resize)> pResize;
private:
	void *m_buf;
	int m_alloc;
	int m_size;
	int m_granul;
};

template<class PTRTYPE> class WDL_TypedBuf {
public:
	PTRTYPE *Get() const { return (PTRTYPE *) m_hb.Get(); }
	int GetSize() const { return m_hb.GetSize()/(unsigned int)sizeof(PTRTYPE); }
    PTRTYPE *Resize(int newsize, bool resizedown = true) { return (PTRTYPE *)m_hb.Resize(newsize*sizeof(PTRTYPE),resizedown); }
private:
	WDL_HeapBuf m_hb;
};


class WDL_Resampler {
public:
	int ResampleOut(WDL_ResampleSample *out, int nsamples_in, int nsamples_out, int nch);
	static MemberFuncPtr<decltype(&WDL_Resampler::ResampleOut)> pResampleOut;
private:
	void BuildLowPass(double filtpos);
	const WDL_SincFilterSample *GetFilterCoeff();
	inline void SincSample2N(WDL_ResampleSample *outptr, const WDL_ResampleSample *inptr, double fracpos, int nch, const WDL_SincFilterSample *filter, int filtsz);
	
	double m_sratein; //WDL_FIXALIGN;
	double m_srateout;
	double m_fracpos;
	double m_ratio;
	double m_filter_ratio;
	float m_filterq, m_filterpos;
	WDL_TypedBuf<WDL_ResampleSample> m_rsinbuf;
	WDL_TypedBuf<WDL_SincFilterSample> m_filter_coeffs;

	void *m_iirfilter;

	int m_filter_coeffs_size;
	int m_last_requested;
	int m_filtlatency;
	int m_samples_in_rsinbuf;
	int m_lp_oversize;

	int m_sincsize;
	int m_filtercnt;
	int m_sincoversize;
	char m_interp;
	char m_feedmode;
	
	
	class Cache {
	public:
		Cache() : coeff(NULL), state(0) {}
		WDL_SincFilterSample *coeff;
		volatile LONG state;
	};
	static Cache m_cache[4][6][6];  // [sincsize][sratein][srateout]
	
	inline static int GetSrateIndex(float srate);
	Cache *GetCache() const;

	
	void _check_struct() {
	#ifdef _WIN64
		static_assert(offsetof(WDL_Resampler, m_feedmode) == 0x89, "struct offset error");
	#else
		static_assert(offsetof(WDL_Resampler, m_feedmode) == 0x75, "struct offset error");
	#endif
	}
};

#endif
