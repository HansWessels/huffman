/*
** testen van een Huffman encoder
** in het bijzonder het package merge algoritme
*/

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h> // https://en.cppreference.com/w/c/types/integer.html#Format_macro_constants
#include <stdlib.h>
#include <string.h>


#define MAX_SYMBOL_SIZE 256
#define MAX_HUFFMAN_LEN 64
#define SYMBOL_SIZE 256

#if MAX_SYMBOL_SIZE <= 256 /* voor max_symbol_size <=256 */
    typedef uint_fast8_t symbol_t;
    typedef int_fast16_t symbol_count_t;
#elif MAX_SYMBOL_SIZE <= 32768 /* voor max_simbol_size <=32768 */
    typedef uint_fast16_t symbol_t;
    typedef int_fast16_t symbol_count_t;
#elif MAX_SYMBOL_SIZE <= 65536 /* voor max_simbol_size <=65536 */
    typedef uint_fast16_t symbol_t;
    typedef int_fast32_t symbol_count_t;
#elif MAX_SYMBOL_SIZE <= (1<<31) /* voor max_simbol_size <=1<<31 */
    typedef uint_fast32_t symbol_t;
    typedef int_fast32_t symbol_count_t;
#elif MAX_SYMBOL_SIZE <= (1<<32) /* voor max_simbol_size <=1<<32 */
    typedef uint_fast32_t symbol_t;
    typedef int_fast64_t symbol_count_t;
#elif 0 /* rest */
    typedef uint_fast64_t symbol_t;
    typedef int_fast64_t symbol_count_t;
#endif


#if MAX_HUFFMAN_LEN <= 8
    typedef uint_fast8_t huffman_t;
#elif MAX_HUFFMAN_LEN <= 16
    typedef uint_fast16_t huffman_t;
#elif MAX_HUFFMAN_LEN <= 32
    typedef uint_fast32_t huffman_t;
#else
    typedef uint64_t huffman_t;
#endif
typedef uint64_t freq_t;


/*
** ARJ CRC32 routines
**
** 2022 Hans Wessels
*/

/*
** this function makes a CRC32 table for ARJ CRC32
** call function with pointer to an unsigned long[256] array
*/

void make_crc32_table(uint32_t crc_table[])
{
  int max_bits=8;
  int bits=1;
  crc_table[0]=0;
  crc_table[1<<(max_bits-1)]=0xEDB88320UL;
  do
  {
    int i=(1<<bits)-1;
    do
    {
      int j=(1<<bits)-1;
      int done=i<<(max_bits-2*bits);
      uint32_t crc=crc_table[done<<bits];
      int offset=(int)((crc&j)<<(max_bits-bits));
      crc>>=bits;
      do
      {
        int tmp=j<<(max_bits-bits);
        crc_table[tmp+done]=crc^crc_table[tmp^offset];
      }
      while(j--);
    }
    while(--i);
    bits+=bits;
  }
  while(bits<max_bits);
}

/*
** this function calculates the ARJ CRC32 checksum over count bytes in data
** call function with count: nuber of bytes to CRC, data: pointer to data bytes
*/

uint32_t crc32(unsigned long count, uint8_t* data, uint32_t crc_table[])
{
	uint32_t crc=(uint32_t)-1;
	while(count--!=0)
	{
		uint8_t c = *data++;
		c ^= (uint8_t) crc;
		crc >>= 8;
		crc ^= crc_table[c];
	}
	crc=~crc;
	return crc;
}

uint32_t crc32init(void)
{
	uint32_t crc=(uint32_t)-1;
	return crc;
}

uint32_t crc32byte(uint8_t c, uint32_t crc_table[], uint32_t crc)
{
	c ^= (uint8_t) crc;
	crc >>= 8;
	crc ^= crc_table[c];
	return crc;
}

int64_t load_file(char* infile, uint8_t** data_in)
{
	FILE* f;
	uint8_t* data;
	int64_t size;
	*data_in=NULL;
	f = fopen(infile, "rb");
	if (f == NULL)
	{
		fprintf(stderr, "File open error %s!\n", infile);
		return -1;
	}
	fseek(f, 0, SEEK_END);
	size = ftell(f);
	data = (uint8_t*)malloc(size+1024);
	if (data == NULL)
	{
		fprintf(stderr, "Malloc error voor file data %s !\n", infile);
		fclose(f);
		return -1;
	}
	fseek(f, 0, SEEK_SET);
	if(fread(data, 1, size, f)!=size)
	{
		fprintf(stderr, "Read error %s\n", infile);
		fclose(f);
		return -1;
	}
	fclose(f);
	*data_in=data;
	return size;
}

void freq_count(const uint8_t *data, int64_t size, freq_t freq[], symbol_count_t symbol_size)
{
    if(symbol_size>(1<<24))
    {
        while(size>=4)
        {
            uint32_t tmp;
            tmp=data[size--];
            tmp<<=8;
            tmp+=data[size--];
            tmp<<=8;
            tmp+=data[size--];
            tmp<<=8;
            tmp+=data[size--];
            freq[tmp]++;
        }
    }
    else if(symbol_size>(1<<16))
    {
        while(size>=3)
        {
            uint32_t tmp;
            tmp=data[size--];
            tmp<<=8;
            tmp+=data[size--];
            tmp<<=8;
            tmp+=data[size--];
            freq[tmp]++;
        }
    }
    else if(symbol_size>(1<<8))
    {
        while(size>=2)
        {
            uint16_t tmp;
            tmp=data[size--];
            tmp<<=8;
            tmp+=data[size--];
            freq[tmp]++;
        }
    }
    else
    {
        while(size>0)
        {
            freq[data[size--]]++;
        }
    }
}

int freq_compare(const void* a_in, const void* b_in)
{
    freq_t* a=(freq_t *)a_in;
    freq_t* b=(freq_t *)b_in;
    if(*a<*b)
    {
        return -1;
    }
    else if(*a>*b)
    {
        return 1;
    }
    return 0;
}

void sort_symbols(symbol_t symbols[MAX_SYMBOL_SIZE], freq_t freq[MAX_SYMBOL_SIZE], symbol_count_t symbol_count)
{
    freq_t pairs[2*MAX_SYMBOL_SIZE];
    symbol_count_t i;
    for(i=0; i<symbol_count; i++)
    {
        pairs[2*i]=freq[symbols[i]];
        pairs[2*i+1]=symbols[i];
    }
    qsort(pairs, symbol_count, 2*sizeof(uint64_t), freq_compare);
    for(i=0; i<symbol_count; i++)
    {
        symbols[i]=(symbol_t)pairs[2*i+1];
    }
}

void walk_tree(symbol_count_t node, symbol_count_t c_left[MAX_SYMBOL_SIZE], symbol_count_t c_right[MAX_SYMBOL_SIZE], int s_len[])
{
    if(node<0)
    {
        walk_tree(c_left[node], c_left, c_right, s_len);
        walk_tree(c_right[node], c_left, c_right, s_len);
    }
    else
    {
        s_len[node]++;
    }
}

void make_huffman_codes(int s_len[], huffman_t huff_codes[], symbol_count_t symbol_count)
{
    int len_count[MAX_HUFFMAN_LEN+1]={0};
    huffman_t huffcode[MAX_HUFFMAN_LEN+1];
    symbol_count_t i;
    i=symbol_count;
    do
    {
        i--;
        len_count[s_len[i]]++;
    } while(i>0);
    huffman_t start_huffcode=0;
    for(i=1; i<=MAX_HUFFMAN_LEN; i++)
    {
        start_huffcode<<=1;
        huffcode[i]=start_huffcode;
        start_huffcode+=len_count[i];
    }
    for(i=0; i<symbol_count; i++)
    {
        if(s_len[i]>0)
        {
            huff_codes[i]=huffcode[s_len[i]];
            huffcode[s_len[i]]++;
        }
        else
        {
            huff_codes[i]=0;
        }
    }
}

int make_huffman_table(int s_len[], huffman_t huff_codes[], const freq_t in_freq[], int max_huff_len, const symbol_count_t symbol_size)
{
    freq_t freq_array[MAX_SYMBOL_SIZE+MAX_SYMBOL_SIZE*MAX_HUFFMAN_LEN]={0};
    symbol_count_t c_left_array[MAX_SYMBOL_SIZE*MAX_HUFFMAN_LEN];
    symbol_count_t c_right_array[MAX_SYMBOL_SIZE*MAX_HUFFMAN_LEN];
    symbol_t symbols[MAX_SYMBOL_SIZE];
    symbol_count_t pairs[MAX_SYMBOL_SIZE];
    freq_t* freq=freq_array+MAX_SYMBOL_SIZE*MAX_HUFFMAN_LEN;
    symbol_count_t* c_left=c_left_array+MAX_SYMBOL_SIZE*MAX_HUFFMAN_LEN;
    symbol_count_t* c_right=c_right_array+MAX_SYMBOL_SIZE*MAX_HUFFMAN_LEN;
    symbol_count_t symbol_count=0;
    symbol_count_t pairs_count;
    symbol_count_t node=-1;
    const uint64_t een=1;
    if(max_huff_len>MAX_HUFFMAN_LEN)
    {
        printf("max_huff_len te groot: %i\n", max_huff_len);
        return -1;
    }
    if(max_huff_len<64)
    {
        if((een<<max_huff_len)<symbol_size)
        {
            printf("(1<<max_huff_len) < symbol_size : 1<<%" PRIX64 " < %" PRIX64 "\n", een<<max_huff_len, symbol_size);
            return -1;
        }
    }
    symbol_count_t i=symbol_size;
    do
    { /* hoeveel symbols zijn er met een freq>0? */
        i--;
        freq[i]=in_freq[i];
        if(freq[i]!=0)
        {
            symbols[symbol_count]=i;
            symbol_count++;
        }
    } while(i>0);
    sort_symbols(symbols, freq, symbol_count);
    pairs_count=symbol_count>>1;
    i=pairs_count;
    do
    { /* eerste ronde, gewoon de gegeven character bij elkaar voegen */
        i--;
        freq[node]=freq[symbols[i*2]]+freq[symbols[i*2+1]];
        c_left[node]=symbols[i*2];
        c_right[node]=symbols[i*2+1];
        pairs[i]=node;
        node--;
    } while(i>0);
    max_huff_len--;
    do
    { /* merge symbols, max_huff_len-1 keer */
        symbol_count_t symbol_pos=symbol_count-1;
        symbol_count_t pair_pos=pairs_count-1;
        freq_t next_symbol_freq=freq[symbols[symbol_pos]];
        freq_t next_pair_freq=freq[pairs[pair_pos]];
        pairs_count=symbol_count+pairs_count;
        if((pairs_count)&1)
        { /* oneven som, waarde met hoogste freq doet niet mee */
            if(next_pair_freq>next_symbol_freq)
            {
                pair_pos--;
                next_pair_freq=freq[pairs[pair_pos]];
            }
            else
            {
                symbol_pos--;
                next_symbol_freq=freq[symbols[symbol_pos]];
            }
        }
        pairs_count>>=1;
        i=pairs_count;
        for(;;)
        { /* maak de nieuwe pairs, exit lus altijd doordat de pairs op zijn */
            freq_t node_freq;
            i--;
            if(next_pair_freq>next_symbol_freq)
            {
                node_freq=next_pair_freq;
                c_left[node]=pairs[pair_pos];
                pair_pos--;
                if(pair_pos<0)
                { /* pairs zijn op, koppel met symbol en exit */
                    node_freq+=next_symbol_freq;
                    freq[node]=node_freq;
                    c_right[node]=symbols[symbol_pos];
                    symbol_pos--;
                    break;
                }
                else
                {
                    next_pair_freq=freq[pairs[pair_pos]];
                }
            }
            else
            {
                node_freq=next_symbol_freq;
                c_left[node]=symbols[symbol_pos];
                symbol_pos--;
                next_symbol_freq=freq[symbols[symbol_pos]];
            }
            if(next_pair_freq>next_symbol_freq)
            {
                node_freq+=next_pair_freq;
                c_right[node]=pairs[pair_pos];
                pair_pos--;
                if(pair_pos<0)
                { /* pairs zijn op, exit */
                    freq[node]=node_freq;
                    break;
                }
                else
                {
                    next_pair_freq=freq[pairs[pair_pos]];
                }
            }
            else
            {
                node_freq+=next_symbol_freq;
                c_right[node]=symbols[symbol_pos];
                symbol_pos--;
                next_symbol_freq=freq[symbols[symbol_pos]];
            }
            freq[node]=node_freq;
            pairs[i]=node;
            node--;
        }
        pairs[i]=node;
        node--;
        do
        { /* koppel de laatste pairs */
            i--;
            freq[node]=freq[symbols[i*2]]+freq[symbols[i*2+1]];
            c_left[node]=symbols[i*2];
            c_right[node]=symbols[i*2+1];
            pairs[i]=node;
            node--;
        } while(i>0);
        max_huff_len--;
    } while(max_huff_len>0);
    memset(s_len, 0, symbol_size*sizeof(s_len[0]));
    i=pairs_count; /* N symbolen levert altidj N-1 pairs op */
    do
    { /* s_len wordt berekend door het aantal keer dat een bepaald symbool in de eeste (symbol_count-1) pairs voorkomt */
        i--;
        walk_tree(pairs[i], c_left, c_right, s_len);
    } while(i>0);
    make_huffman_codes(s_len, huff_codes, symbol_size);
    return 0;
}

uint64_t encode(const uint8_t *data, const int64_t size, const symbol_count_t symbol_size, const int max_huff_len)
{
    freq_t freq[SYMBOL_SIZE]={0};
    int s_len[SYMBOL_SIZE];
    huffman_t huffcodes[SYMBOL_SIZE];
    freq_count(data, size, freq, symbol_size);
    make_huffman_table(s_len, huffcodes, freq, max_huff_len, symbol_size);
    {
        uint64_t totaal=0;
        symbol_count_t i;
        for(i=0; i<symbol_size; i++)
        {
            totaal+=freq[i]*s_len[i];
        }
        return totaal;
    }
}

int main(int argc, char* argv[])
{
	uint32_t crc_table[256];
	int i;
	i=1;
	make_crc32_table(crc_table);
	uint64_t totaal=0;
	while(i<argc)
	{
		uint32_t crc;
		int64_t size;
		uint8_t* data;
		const int symbol_size=256;
		const int max_huff_len=15;
		printf("Loading file %s\n",argv[i]);
		size=load_file(argv[i], &data);
		if(data==NULL)
		{
			printf("Error loading %s\n", argv[i]);
			return -1;
		}
		totaal+=encode(data, size, symbol_size, max_huff_len);
		free(data);
		i++;
	}
	printf("Totaal = %lu\n", totaal);
    return 0;
}
