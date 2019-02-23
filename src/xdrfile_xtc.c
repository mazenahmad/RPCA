/* -*- mode: c; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- 
 *
 * $Id$
 *
 * Copyright (c) 2009-2014, Erik Lindahl & David van der Spoel
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
 
#include <stdlib.h>
#include <math.h>
#include "xdrfile.h"
#include "xdrfile_xtc.h"
#include "arrays.h"
	
#define MAGIC 1995
	char *xdr_message[exdrNR] = {
		"OK",
		"Header",
		"String",
		"Double",
		"Integer",
		"Float",
		"Unsigned integer",
		"Compressed 3D coordinate",
		"Closing file",
		"Magic number",
		"Not enough memory",
		"End of file",
		"File not found"
	};



static int xtc_header(XDRFILE *xd,int *natoms,int *step,float *time,int bRead)
{
	int result,magic,n=1;
	
	/* Note: read is same as write. He he he */
	magic  = MAGIC;
	if ((result = xdrfile_write_int(&magic,n,xd)) != n)
		{
			if (bRead)
				return exdrENDOFFILE;
			else
				return exdrINT;
		}
	if (magic != MAGIC)
		return exdrMAGIC;
	if ((result = xdrfile_write_int(natoms,n,xd)) != n)
		return exdrINT;
	if ((result = xdrfile_write_int(step,n,xd)) != n)
		return exdrINT;
	if ((result = xdrfile_write_float(time,n,xd)) != n)
		return exdrFLOAT;
	
	return exdrOK;
}

static int xtc_coord(XDRFILE *xd,int *natoms,float box[3][3],rvec *x,float *prec, int bRead)
{
	int i,j,result;
    
	/* box */
	result = xdrfile_read_float(box[0],DIM*DIM,xd);
	if (DIM*DIM != result)
		return exdrFLOAT;
	else 
		{
			if (bRead)
				{
					result = xdrfile_decompress_coord_float(x[0],natoms,prec,xd); 
					if (result != *natoms)
						return exdr3DX;
				}
			else
				{
					result = xdrfile_compress_coord_float(x[0],*natoms,*prec,xd); 
					if (result != *natoms)
						return exdr3DX;
				}
		}
	return exdrOK;
}


static int read_xtc( XDRFILE *xd, int *natoms,int *step,float *time, matrix box,rvec *x)
/* Read subsequent frames */
{
	int result;
	float prec;
	if ((result = xtc_header(xd,natoms,step,time,TRUE)) != exdrOK)
		return result;
	if ((result = xtc_coord(xd,natoms,box,x,&prec,1)) != exdrOK)
		return result;

	return exdrOK;
}



int write_xtc(XDRFILE *xd,int natoms,int step,float time, matrix box,rvec *x)
/* Write a frame to xtc file */
{
	int result;
	float prec=1000;
	if ((result = xtc_header(xd,&natoms,&step,&time,FALSE)) != exdrOK)
		return result;

	if ((result = xtc_coord(xd,&natoms,box,x,&prec,0)) != exdrOK)
		return result;
  
	return exdrOK;
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

int read_xtc_natoms(const char *fn,int *natoms)
{
	XDRFILE *xd;
	int step,result;
	float time;

	xd = xdrfile_open(fn,"r");
	if (NULL == xd)
	{
		fprintf(stderr,"File %s can not be opened\n", fn);
		return 0;
	}
	result = xtc_header(xd,natoms,&step,&time,TRUE);
	xdrfile_close(xd);
	if(result == exdrOK )
	{
		return 1;
	}
		else
	{
		fprintf(stderr, "%s\n", xdr_message[result]);
		return 0;
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

int Read_first_xtc( const char *fn,xtc_t *f, float start)
{
	int result;

	read_xtc_natoms(fn,&(f->natoms));
	f->xd = xdrfile_open(fn,"r");
	if (NULL == f->xd)
	{
		PRINT_ERR_OPENFILE( fn);
	}
	ALLO_ONE_DIM_ARRAY(f->x, rvec,(f->natoms) );

	do
	{
		result = read_xtc(f->xd,&(f->natoms),&(f->step),&(f->time),f->box,f->x);
	} while((result == exdrOK ) && ( f->time) < start );
	fprintf(stderr,"Start reading the traj. at %.3f ps\n",f->time );

	if(result == exdrOK )
	{
		fprintf(stderr,"Reading frame at time: ");
		return 1;
	}
	else if (result == exdrENDOFFILE )
	{
		xdrfile_close(f->xd);
		fprintf(stderr, "No more frames to read after the start\n");
		return 0;
	}
	else
	{
		fprintf(stderr, "%s\n", xdr_message[result]);
		return 0;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

int Read_next_xtc( xtc_t *f, float end)
{
	int result = read_xtc(f->xd,&(f->natoms),&(f->step),&(f->time),f->box,f->x);
	if(end <= 0) end = f->time + 100;
	if(result == exdrOK  && (f->time) < end  )
	{
		if( f->time > 1000 )
		{
			if(fmod(f->time ,250 )< 0.01 )
			{
				//fprintf(stderr,"%14.3f ps",f->time );
				fprintf(stderr,"%14.3f ns",(f->time)/1000 );
				fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
			}
		}
		else
		{
			fprintf(stderr,"%14.3f ps",f->time );
			fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
		}
		return 1;
	}
	else if ( (f->time) >= end || (result == exdrENDOFFILE) )
	{
		xdrfile_close(f->xd);
		fprintf(stderr,"\nLast read frame at %.1f ps\n",f->time );
		return 0;
	}
	else
	{
		fprintf(stderr, "%s\n", xdr_message[result]);
		return 0;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////
void reset_frame( xtc_t *f)
{
	f->natoms =0;
	if (f->x ) free(f->x);
	f->step=0;
	f->time =-1;
	int i,d;
	for(i=0;i<3 ;i++) for(d=0;d<3 ;d++) f->box[i][d]=0;
	//if(f->xd) xdrfile_close(f->xd);
	f->xd= NULL;
}
////////////////////////////////////////////////////////////////////////////////////////////////
