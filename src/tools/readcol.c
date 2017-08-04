/* ------------------------------------------------------------------ */
/* WRITTEN BY: Mark den Brok                                          *
 * COMMENTS TO: denbrok@astro.rug.nl                                  *
 *                                                                    *
 * LAST UPDATE: 22-05-2011                                            *
 *                                                                    *
 *          THIS MODULE IS DANGEROUS! DON'T USE IT!                   *
 *        unless of course you know what you 're doing...             *
 *        (the problem is that va_list is platform dependent          *
 *     in x86 gcc is just a stack pointer (which is easy to change)   *
 *         however, on AMD gcc passes it by ref. (Actually            *
 *       in AMD the first n params are spread out over registers))    *
 *                                                                    *
 * This is my own version of IDL's readcol, but different syntax,     *
 * first give file, char* with flags, format, and then all pointers   *
 * Writing this in x86 assembly is a piece of cake, however, in       *
 * c this is a real mess - but as I don't want to add any asm here    *
 * I'm keeping it for now like this.                                  *
 *                                                                    *
 * AND ANOTHER WARNING: NOT ALL FORMATS ARE iMPLEMENTED.              */
/* ------------------------------------------------------------------ */

/*
This code is released under a BSD 2-clause license.
If you use this code for your research, please cite:
Watkins et al. 2013, MNRAS, 436, 2598
"Discrete dynamical models of omega Centauri"
http://adsabs.harvard.edu/abs/2013MNRAS.436.2598W
*/



#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>


void readcol_add_to_end( int, int, char *, va_list );
void sort_stuff_readcol( int, char *, va_list );
int nargs_readcol( char * );

int readcol( char *file, char *flags, char *fmt, ... ) {
    
    char *q;
    int i, quiet, len, ntotal;
    char skipsym;       //skip line starting with ... 
    int bufsize;        //max length of row
    int filelen = 10;   //number of rows in file
    char *buffer;
    va_list ap, ap2;
    FILE *fp;
    
    
    // ---------------------------------
    
    
    // read input parameter string
    
    bufsize = 100;
    q = flags; //for historic reasons, this q is here. for similar reasons
    // I'm using a string for passing params instead of a struct.
    
    // parameter string length
    len = 0;
    for ( i = 0; q[i] != '\0'; i++ ) len++;
    
    skipsym = '\0';
    quiet = 1;
    for ( i = 0; i < ( len - 3 ); i++ ) {
        
        //test for skipsym
        if ( i < ( len - 9 ) ) {
            if ( strncmp( q + i, "skipsym=", 8 ) == 0 ) {
                i += 8;
                skipsym = q[i];
            }
        }
        
        //test for verbose
        if ( i < ( len - 5 ) ) {
            if ( strncmp( q + i, "quiet", 5 ) == 0 ) {
                i += 4;
                quiet = 1;
            }
        }
        
        //test for bufsize
        if ( i < ( len - 9 ) ) {
            if ( strncmp( q + i, "bufsize=", 8 ) == 0 ) {
                i += 8;
                sscanf( q + i, "%i", &bufsize );
            }
        }
        
        //test for filelen
        if ( i < ( len - 9 ) ) {
            if ( strncmp( q + i, "filelen=", 8 ) == 0 ) {
                i += 8;
                sscanf( q + i, "%i", &filelen );
            }
        }
        
    }
    
    if ( !quiet ) {
        printf( "file: %s\n", file );
        printf( "bufsize: %i\n", bufsize );
        printf( "filelen: %i\n", filelen );
        if ( skipsym=='\0' ) printf( "skipsym: '\\0'\n" );
        else printf( "skipsym: '%c'\n", skipsym );
    }
    
    buffer = (char *) malloc( ( bufsize + 1 ) * sizeof( char ) );
    
    
    // ---------------------------------
    
    
    //  HERE STARTS THE REAL WORK
    
    ntotal = 0;
    fp = fopen( file, "r" );
    while( fgets( buffer, bufsize, fp ) != NULL ) {
        if ( buffer[0] != skipsym && buffer[0] != '\n' ) {
            
            va_start( ap, fmt );
            va_copy( ap2, ap );
            va_end( ap );
            
            va_start( ap, fmt );
            if ( vsscanf( buffer, fmt, ap ) > 0 ) {
                readcol_add_to_end( ntotal, filelen, fmt, ap2 );
                ntotal++;
            }
            va_end( ap );
            
        }
    }
    fclose( fp );
    
    va_start( ap, fmt );
    //Keats was right:  A thing of beauty is a joy forever!
    sort_stuff_readcol( filelen, fmt, ap );
    va_end( ap );
    
    if ( !quiet ) {
        printf( "%i valid lines read\n", ntotal );
    }
    
    free( buffer );
    
    return ntotal;
}




void readcol_add_to_end( int ntotal, int filelen, char *fmt, va_list ap ) {
    /*Programming this in c is diff because parameters are protected for
        editing. Instead of changing the pointers in the parameter list,
        the pointers stay constant, but here I copy the content of ..[0] to
        the 'right' location in the array. Note that it will give an array
        which is upside down. We have to correct for this in the end. */
    
    const char *fp;
    char **s, **st;
    char *c, *ct;
    int *i, *it;
    double *d, *dt;
    float *f, *ft;
    
    for ( fp = fmt; *fp != '\0'; fp++ ) {
        
        if ( *fp == '%' )
            switch ( *++fp ) {
                
                case 'l':
                    switch( *++fp ) {
                        
                        case 'f':
                            d = va_arg( ap, double * );
                            dt = d + ( filelen - ntotal - 1 );
                            *dt = *d;
                            break;
                        
                    }
                    break;
                
                case 's':
                    s = va_arg( ap, char ** );
                    st = s + ( filelen - ntotal - 1 );
                    *st = *s;
                    break;
                
                case 'c':
                    c = va_arg( ap, char * );
                    ct = c + ( filelen - ntotal - 1 );
                    *ct = *c;
                    break;
                
                case 'i':
                    i = va_arg( ap, int * );
                    it = i + ( filelen - ntotal - 1 );
                    *it = *i;
                    break;
                
                case 'f':
                    f = ( float * ) va_arg( ap, float * );
                    ft = f + ( filelen - ntotal - 1 );
                    *ft = *f;
                    break;
                
            }
        
    }
    
    return;
}


void sort_stuff_readcol( int filelen, char *fmt,va_list ap ) {
    /* Since all read values are stored from the back of the array to the
        front, we have to reverse them again. */
    
    const char *fp;
    char **s, **st, *sq;
    char *c, *ct, cq;
    int *i, *it, iq;
    double *d, *dt, dq;
    float *f, *ft, fq;
    int j, jmax;
    jmax = filelen / 2;
    
    for ( fp = fmt; *fp != '\0'; fp++ ) {
        if ( *fp == '%' )
            switch ( *++fp ) {
                
                case 'l':
                    switch ( *++fp ) {
                        
                        case 'f':
                            d = va_arg( ap, double * );
                            dt = d + ( filelen - 1 );
                            for ( j = 0; j < jmax; j++ ) {
                                dq = *d;
                                *d = *dt;
                                *dt = dq;
                                d++;
                                dt--;
                            }
                            break;
                        
                    }
                    break;
                    
                case 's':
                    s = va_arg( ap, char ** );
                    st = s + ( filelen - 1 );
                    for ( j = 0; j < jmax; j++ ) {
                        sq = *s;
                        *s = *st;
                        *st = sq;
                        s++;
                        st--;
                    }
                    break;
                
                case 'c':
                    c = va_arg( ap, char * );
                    ct = c + ( filelen - 1 );
                    for ( j = 0; j < jmax; j++ ) {
                        cq = *c;
                        *c = *ct;
                        *ct = cq;
                        c++;
                        ct--;
                    }
                    break;
                
                case 'i':
                    i = va_arg( ap, int * );
                    it = i + ( filelen - 1 );
                    for ( j = 0; j < jmax; j++ ) {
                        iq = *i;
                        *i = *it;
                        *it = iq;
                        i++;
                        it--;
                    }
                    break;
                
                case 'f':
                    f = va_arg( ap, float * );
                    ft = f + ( filelen - 1 );
                    for ( j = 0; j < jmax; j++ ) {
                        fq = *f;
                        *f = *ft;
                        *ft = fq;
                        f++;
                        ft--;
                    }
                    break;
                
            }
        
    }
    
    return;
    
}



int nargs_readcol( char *fmt ) {
    
    const char *fp;
    
    int j;
    j = 0;
    
    for ( fp = fmt; *fp != '\0'; fp++ )
        if ( *fp == '%' && ( *( fp + 1 ) != '*' ) )
            j++;
    
    return j;
    
}
