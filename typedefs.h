#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#ifdef _WIN32
typedef	unsigned short	Uint16;
typedef	short			int16;
typedef	unsigned int	Uint32;
typedef int				int32;
typedef double			float32;
typedef double			float64;
#else
typedef	unsigned int	Uint16;
typedef	int				int16;
typedef	unsigned long	Uint32;
typedef long			int32;
typedef float			float32;
typedef long double		float64;
#endif

#endif//_TYPEDEFS_H_