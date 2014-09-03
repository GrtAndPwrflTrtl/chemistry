/****************************************************************************
*
* Title:     bstream.c
*
* Synopsis:  This module contains routines to perform the following
*            processing:
*
* 1.  I/O for bit streams
*
*     NOTE (warning): implementation is done with an unsigned char type for
*     the stream.  Therefore, if any manipulation of the stream is done
*     independent of the functions provided here, endian problems may occur.
*
*     Optimization (gcc -03 inlining, etc) test showed a 26% decrease in
*     execution time (150000 -> 110000) with a slight increase in code size
*     (30680 -> 31200 bytes) : 1.67%
*
* Created by: C. Alvin
*
****************************************************************************/

/*
 * Include Files : std C libs
 */

#include <stdio.h>
#include <stdlib.h> /* malloc */
#include <string.h> /* memset */


/*
 * Include Files : project related
 */

#include "bstream.h"


/***************** STREAM INDEPENDENT FUNCTIONS ************************/

/*#######################################################################
 *
 *      Synopsis: resets a bit stream so that it can be used as input
 *                from the beginning of the stream.  Should only be used
 *                if the intent is input from a stream. A new stream
 *                must be created every time for safety.
 *                Function NOT meant as a clearing mechanism.
 *
 *      Parameters: a bit stream
 *
 *      Return Value: NONE
 *
 *#######################################################################*/

void resetBitStream(bstream* out)
{
  out->current = out->stream;
  out->mask = 0x80;
  out->end = out->stream + out->stream_length;
}


/*#######################################################################
 *
 *      Synopsis: a C destructor for a bit stream
 *
 *      Parameters: a bit stream
 *
 *      Return Value: NONE
 *
 *#######################################################################*/

void freeBitStream(bstream* s)
{
  if (s->stream) free(s->stream);
  s->end = s->current = s->stream = NULL;
}

/*#######################################################################
 *
 *      Synopsis: clears and resets a bit stream so that it may be
 *
 *      Parameters: a bit stream
 *
 *      Return Value: NONE
 *
 *#######################################################################*/

void clearBitStream(bstream* s)
{
  memset(s->stream, 0, s->MAX_SIZE); /* clear the stream to be all Os */
  s->end = s->current = s->stream; /* s->end useless until reset for reading */
  s->mask = 0x80;
  s->stream_length = 0;
  s->garbage_bits = 0;
  /* s->MAX_SIZE  remains unmodified */
}

/*******************  BIT  STREAM  FUNCTIONS   *********************/

/*##########################################################################
 *
 *      Synopsis: a C constructor for a bit stream
 *
 *      Parameters: the length of the desired bit stream (in bytes)
 *
 *      Return Value: a newly allocated and initialized bit stream
 *
 *##########################################################################*/

bstream* openBitStream(unsigned length)
{
  bstream* out = (bstream*)malloc(sizeof(bstream));
  out->current = out->stream =
             (unsigned char*)malloc(sizeof(unsigned char) * length + 1);

  /* zero out stream */
  memset(out->stream, 0, sizeof(unsigned char) * length + 1);

  out->end = out->current; /* useless until reset for reading */
  out->mask = 0x80;
  out->stream_length = 0;
  out->garbage_bits = 0;
  out->MAX_SIZE = length;

  return out;
}

/*##########################################################################
 *
 *      Synopsis: adds one bit to the stream
 *
 *      Restrictions: requires that the bit to write is found in the least
 *                     significant bit position (0x0 or 0x1) (big endian)
 *
 *      Note: want to inline this function if possible
 *
 *      Parameters: a bit stream and bit to write
 *
 *      Return Value: NONE
 *
 *##########################################################################*/

void putBit(bstream* out, unsigned char bit)
{
  if (!out->stream)
  {
    fprintf(stderr, "ERROR: stream has not been initialized.\n");
    return;
  }

  /* if a 1, mask in 1 (at current bit stream mask position)  */
  if (bit & 0x1)
  {
    *out->current |= out->mask;
  }

  out->mask >>= 1;

  /* Check if the byte is filled and can be written */
  if (!out->mask /* == 0 */)
  {
    if (out->stream_length == (int)out->MAX_SIZE)
    {
      fprintf(stderr, "putBits::stream size exceeded\n");
    }

    /* move ahead work byte and reset mask */
    *++out->current = 0;
    out->mask = 0x80;

    /* increment stream length */
    out->stream_length++;
  }
}

/*##########################################################################
 *
 *      Synopsis: adds a specified number of bits to the stream.
 *                Does not call putBit in attempt to execute faster.
 *
 *      Restrictions: requires that the bits to write is found in the least
 *                     significant bit position (0x0 or 0x1) (big endian)
 *
 *      Note: want to inline this function if possible
 *
 *      Parameters: a bit stream, value to write, and # bits to write
 *
 *      Return Value: NONE
 *
 *##########################################################################*/

void putBits(bstream* out, unsigned bits, unsigned short num_bits)
{
  if (num_bits == 0)
  {
    return;
  }

  if (num_bits == 1)
  {
    putBit(out, bits);
    return;
  }

  unsigned mask = 1 << (num_bits - 1);
  unsigned short bit_counter = num_bits;

  if (!out->stream)
  {
    fprintf(stderr, "ERROR: stream has not been initialized.\n");
    return;
  }

  /* error handle */
  if (bit_counter > 32)
  {
    fprintf(stderr, "putBits::attempting to write more than 32 bits");

    /* reset values appropriately */
    bit_counter = 32;
    mask = 0x80000000;
  }

  /* finish filling the current byte */
  while (out->mask != 0x80 && bit_counter)
  {
    putBit(out, bits & mask ? 1 : 0);
    mask >>= 1;
    bit_counter--;
  }

  /* if the mask is at the beginning of a byte and bits have not been read */
  /* in this function, then advance stream pointer */
  /* this is performed in putBit */
  /* if (out->mask == 0x80 && bit_counter != num_bits) out->current++; */

  /* write byte by byte : speed things up */
  while (bit_counter >= 8)
  {
    if (out->stream_length == (int)out->MAX_SIZE)
    {
      fprintf(stderr, "putBits:: stream size exceeded\n");
    }

    /* mask out correct 8 bits */
    *out->current++ = (unsigned char)
                      (((bits & mask)  | (bits & (mask - 1))) >> (bit_counter - 8));
    out->stream_length++;
    bit_counter -= 8;
    mask >>= 8;
  }

  /* write remaining bits to the current byte */
  while (mask)
  {
    if (bits & mask) *out->current |= out->mask;
    out->mask >>= 1;
    mask >>= 1;
  }
}

/*##########################################################################
 *
 *      Synopsis: closes the bit stream and cleans up the remaining byte by
 *                adding garbage bits.  Also, resets the stream so it can
 *                immediately be read from.
 *
 *      Parameters: a bit stream
 *
 *      Return Value: the number of garbage bits added to the end of the stream
 *
 *##########################################################################*/
unsigned short closeBitStream(bstream* out)
{
  /* if there are leftover bits, output them to stream */
  /* also pad and track garbage bits added */
  while (out->mask != 0x80)
  {
    putBit(out, 0); /* putBit updates out->mask */
    out->garbage_bits++;
  }

  if (out->garbage_bits > 7)
  {
    fprintf(stderr,
            "closeBitStream::peculiar number of garbage bits (> 7)");
  }

  /* reset so ok to read from */
  resetBitStream(out) ;

  return out->garbage_bits;
}

/*##########################################################################
 *
 *      Synopsis: retrieves a specified number of bits from the stream
 *                leaving the stream unmodified.
 *
 *      Restrictions: Only 32 bits or less are allowed to be peeked at.
 *
 *      Parameters: a bit stream, # bits to read
 *
 *      Return Value: an unsigned value containing the desired number of bits
 *
 *##########################################################################*/

unsigned peekBits(bstream* in, unsigned short num_bits)
{
  unsigned char* curr = in->current;
  unsigned char msk = in->mask;
  unsigned value = 0;
  unsigned short i;

  if (curr == in->end)
  {
    fprintf(stderr, "ERROR: peeking at end of stream.\n");
    return 0;
  }

  /* mask in number of desired bits - 1 ; last bit is special */
  for (i = 0; i < num_bits - 1; i++)
  {
    if (*curr & msk) value |= 0x1;
    value <<= 1;
    msk >>= 1;

		/* check if need to go to next byte */
    if (!msk)
    {
      msk = 0x80;
      curr++;
      if (curr == in->end) return value; 
    }
  }

  /* handle last (possibly 32nd) bit */
  if (num_bits)
  {
    if (*curr & msk) value |= 0x1;
  }

  return value;
}

Peek from fron

/*####################################################################^
 *
 *      Synopsis: retrieves a bit from the stream.
 *
 *      Parameters: a bit stream
 *
 *      Return Value: an unsigned char with the specified bit in the least
 *                    significant position.
 *
 *##########################################################################*/

unsigned char getBit(bstream* in)
{
  unsigned val;

  if (in->current == in->end)
  {
    fprintf(stderr, "getBit::accessing end of input bit stream\n");
    return -1; /* OxFF */
  }

  if ((in->current + 1) == in->end && (in->mask >> in->garbage_bits) == 0)
  {
    fprintf(stderr, "getBit::accessing garbage bits in bit stream\n");
  }

  val = *in->current & in->mask;

  /* shift mask and check if should be reset and current pointer incremented */
  if ((in->mask >>= 1) == 0)
  {
    in->mask = 0x80;
    in->current++;
  }

  return val ? 0x1 : 0x0;
}

/*##########################################################################
 *
 *      Synopsis: retrieves a specified number of bits from the stream.
 *                Calls getBit() for each, can be made faster by removing this.
 *
 *      Restrictions: Only 32 bits or less are allowed to be read
 *                       ??? Endian issues ???
 *
 *      Parameters: a bit stream, # bits to read
 *
 *      Return Value: an unsigned value containing the desired number of bits
 *
 *##########################################################################*/

unsigned getBits(bstream* in, unsigned short num_bits)
{
  unsigned value = 0; /* return value */
  unsigned mask = 1 << (num_bits - 1); /* return value mask */
  unsigned char c;    /* dummy to receive from getBit () */
  unsigned short bit_counter = num_bits;

  if (in->current == in->end)
  {
    fprintf(stderr, "getBits::end of input bit stream");
    return  -1; /* OxFFFFFFFF */
  }

  if (bit_counter > 32)
  {
    fprintf(stderr, "getBits::attempting to read more than 32 bits");

    bit_counter = 32;
    mask = 1 << (bit_counter - 1);
  }

  /* read rest of current byte */
  while (in->mask != 0x80 && bit_counter)
  {
    if ((c = getBit(in)) == (unsigned char) -1) /* end of stream */
    { 
      return value & 0xffffffff;
    }

    /* or in, shift, and update bit_counter */
    value |= c ? mask : 0;
    mask >>= 1 ;
    bit_counter--;
  }

  /* if the mask is at the beginning of a byte and bits have not been read */
  /* in this function,  then advance stream pointer */
  /* if (in->mask == 0x80 && bit_counter != num_bits) in->current++; */

  /* read byte by byte : speed things up */
  while (bit_counter >= 8)
  {
    if (in->current == in->end)
    {
      fprintf(stderr, "getBits::end of input bit stream");
      return -1; /* OxFFFFFFFF */
    }

    /* mask in correct 8 bits */
    value |= (unsigned)(*in->current++) << (bit_counter - 8);
    bit_counter -= 8;
    mask >>= 8;
  }

  /* read remaining bits from current byte */
  while (mask)
  {
    if (*in->current & in->mask)
    {
      value |= mask;
    }
    in->mask >>= 1;
    mask >>= 1;
  }

  return value;
}

/* NOTE for the following functions: */
/* To insure that these functions are inlined and thus reducing execution */
/* time, make sure that the gcc -03 flag is turned on during compilation */

/*##########################################################################
 *
 *     Synopsis: Prints the contents of the stream in hex. 
 *
 *     Parameters: a bit stream
 *
 *     Return Value: N/A
 *
 *##########################################################################*/

void dumpStream(bstream* stream)
{
  int i;

  for (i = 0; i < stream->stream_length; i++)
  {
    printf("%08d\n", stream->stream[i]);
  }
}
