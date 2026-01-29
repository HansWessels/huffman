/*
** testen van een Huffman encoder
** in het bijzonder het package merge algoritme
*/

#include <stdio.h>
#include <stdint.h>
#include <inttypes.h> // https://en.cppreference.com/w/c/types/integer.html#Format_macro_constants
#include <stdlib.h>
#include <string.h>

uint64_t rdtsc(void);
#if 01
uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}
#endif


#define MAX_SYMBOL_SIZE 512
#define MAX_HUFFMAN_LEN 64
#define SYMBOL_SIZE 256

#if MAX_SYMBOL_SIZE <= 256 /* voor max_symbol_size <=256 */
    typedef uint8_t symbol_t;
    typedef int16_t symbol_count_t;
#elif MAX_SYMBOL_SIZE <= 32768 /* voor max_simbol_size <=32768 */
    typedef uint16_t symbol_t;
    typedef int16_t symbol_count_t;
#elif MAX_SYMBOL_SIZE <= 65536 /* voor max_simbol_size <=65536 */
    typedef uint16_t symbol_t;
    typedef int32_t symbol_count_t;
#elif MAX_SYMBOL_SIZE <= (1<<31) /* voor max_simbol_size <=1<<31 */
    typedef uint32_t symbol_t;
    typedef int32_t symbol_count_t;
#elif MAX_SYMBOL_SIZE <= (1<<32) /* voor max_simbol_size <=1<<32 */
    typedef uint32_t symbol_t;
    typedef int64_t symbol_count_t;
#elif 0 /* rest */
    typedef uint64_t symbol_t;
    typedef int64_t symbol_count_t;
#endif


#if MAX_HUFFMAN_LEN <= 8
    typedef uint8_t huffman_t;
#elif MAX_HUFFMAN_LEN <= 16
    typedef uint16_t huffman_t;
#elif MAX_HUFFMAN_LEN <= 32
    typedef uintt32_t huffman_t;
#else
    typedef uint64_t huffman_t;
#endif
typedef uint64_t freq_t;

uint64_t time_huffman=0;
uint64_t time_sort=0;
uint64_t time_walktree=0;


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
    symbol_count_t i;
    for(i=0; i<symbol_size; i++)
    {
        freq[i]=0;
    }
    if(symbol_size>(1<<24))
    {
        while(size>=4)
        {
            uint32_t tmp;
            tmp=data[--size];
            tmp<<=8;
            tmp+=data[--size];
            tmp<<=8;
            tmp+=data[--size];
            tmp<<=8;
            tmp+=data[--size];
            freq[tmp]++;
        }
    }
    else if(symbol_size>(1<<16))
    {
        while(size>=3)
        {
            uint32_t tmp;
            tmp=data[--size];
            tmp<<=8;
            tmp+=data[--size];
            tmp<<=8;
            tmp+=data[--size];
            freq[tmp]++;
        }
    }
    else if(symbol_size>(1<<8))
    {
        while(size>=2)
        {
            uint16_t tmp;
            tmp=data[--size];
            tmp<<=8;
            tmp+=data[--size];
            freq[tmp]++;
        }
    }
    else
    {
        while(size>0)
        {
            freq[data[--size]]++;
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
        pairs[2*i]=freq[i];
        pairs[2*i+1]=symbols[i];
    }
    qsort(pairs, symbol_count, 2*sizeof(uint64_t), freq_compare);
    for(i=0; i<symbol_count; i++)
    {
        freq[i]=pairs[2*i];
        symbols[i]=(symbol_t)pairs[2*i+1];
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
    freq_t freq_array[MAX_SYMBOL_SIZE+1]={0};
    symbol_count_t tree_array[MAX_HUFFMAN_LEN*MAX_SYMBOL_SIZE*2];
    symbol_count_t* tree=tree_array;
    symbol_t symbols[MAX_SYMBOL_SIZE];
    freq_t pairs_freq[MAX_SYMBOL_SIZE];
    freq_t* freq=freq_array+1;
    symbol_count_t symbol_count=0;
    symbol_count_t pairs_count;
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
        if(in_freq[i]!=0)
        {
            freq[symbol_count]=in_freq[i];
            symbols[symbol_count]=i;
            symbol_count++;
        }
    } while(i>0);
    if(symbol_count<3)
    {
        return -1;
    }
    {
        uint64_t time_start=rdtsc();
        sort_symbols(symbols, freq, symbol_count);
        time_sort+=rdtsc()-time_start;
    }
    pairs_count=symbol_count>>1;
    i=pairs_count;
    do
    { /* eerste ronde, gewoon de gegeven character bij elkaar voegen */
        i--;
        pairs_freq[i]=freq[i*2]+freq[i*2+1];
        tree[i*2]=i*2;
        tree[i*2+1]=i*2+1;
    } while(i>0);
    max_huff_len--;
    do
    { /* merge symbols, max_huff_len-1 keer */
        symbol_count_t symbol_pos=symbol_count-1;
        symbol_count_t pair_pos=pairs_count-1;
        freq_t next_symbol_freq=freq[symbol_pos];
        freq_t next_pair_freq=pairs_freq[pair_pos];
        tree+=symbol_count*2;
        pairs_count=symbol_count+pairs_count;
        if((pairs_count)&1)
        { /* oneven som, waarde met hoogste freq doet niet mee */
            if(next_pair_freq>next_symbol_freq)
            {
                pair_pos--;
                next_pair_freq=pairs_freq[pair_pos];
            }
            else
            {
                symbol_pos--;
                next_symbol_freq=freq[symbol_pos];
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
                tree[i*2+1]=-1;
                pair_pos--;
                if(pair_pos<0)
                { /* pairs zijn op, koppel met symbol en exit */
                    node_freq+=next_symbol_freq;
                    pairs_freq[i]=node_freq;
                    tree[i*2]=symbol_pos;
                    symbol_pos--;
                    break;
                }
                else
                {
                    next_pair_freq=pairs_freq[pair_pos];
                }
            }
            else
            {
                node_freq=next_symbol_freq;
                tree[i*2+1]=symbol_pos;
                symbol_pos--;
                next_symbol_freq=freq[symbol_pos];
            }
            if(next_pair_freq>next_symbol_freq)
            {
                node_freq+=next_pair_freq;
                tree[i*2]=-1;
                pair_pos--;
                if(pair_pos<0)
                { /* pairs zijn op, exit */
                    pairs_freq[i]=node_freq;
                    break;
                }
                else
                {
                    next_pair_freq=pairs_freq[pair_pos];
                }
            }
            else
            {
                node_freq+=next_symbol_freq;
                tree[i*2]=symbol_pos;
                symbol_pos--;
                next_symbol_freq=freq[symbol_pos];
            }
            pairs_freq[i]=node_freq;
        }
        do
        { /* koppel de laatste pairs */
            i--;
            pairs_freq[i]=freq[i*2]+freq[i*2+1];
            tree[i*2]=i*2;
            tree[i*2+1]=i*2+1;
        } while(i>0);
        max_huff_len--;
    } while(max_huff_len>0);
    memset(s_len, 0, symbol_size*sizeof(s_len[0]));
    {
        uint64_t time_start=rdtsc();
        symbol_count_t i;
        for(i=-1; i<symbol_count; i++)
        {
            freq[i]=0;
        }
        i=2*(symbol_count-1);  /* N symbolen levert altidj N-1 pairs op */
        do
        { /* s_len wordt berekend door het aantal keer dat een bepaald symbool in de eeste (symbol_count-1) pairs voorkomt */
            freq[-1]=0;
            do
            {
                i--;
                freq[tree[i]]++;
            } while(i>0);
            i=2*freq[-1];
            tree-=symbol_count*2;
        } while(i>0);
        for(i=0; i<symbol_count; i++)
        {
            s_len[symbols[i]]=(int)freq[i];
        }
        time_walktree+=rdtsc()-time_start;
    }
    make_huffman_codes(s_len, huff_codes, symbol_size);
    return 0;
}

uint64_t encode(const uint8_t *data, const int64_t size, const symbol_count_t symbol_size, const int max_huff_len, const uint64_t block_size)
{
    freq_t freq[SYMBOL_SIZE]={0};
    int s_len[SYMBOL_SIZE];
    huffman_t huffcodes[SYMBOL_SIZE];
    uint64_t pos=0;
    uint64_t delta=block_size;
    uint64_t totaal=0;
    do
    {
        if((pos+delta)>size)
        {
            delta=size-pos;
        }
        if(delta>1)
        {
            freq_count(data+pos, delta, freq, symbol_size);
            uint64_t time_start=rdtsc();
            make_huffman_table(s_len, huffcodes, freq, max_huff_len, symbol_size);
            time_huffman+=rdtsc()-time_start;
            {
                symbol_count_t i;
                for(i=0; i<symbol_size; i++)
                {
                    totaal+=freq[i]*s_len[i];
                }
                /* sanity check Huffcodes */
                {
                    symbol_count_t i;
                    huffman_t max_code=0;
                    for(i=0; i<symbol_size; i++)
                    {
                        if(huffcodes[i]>max_code)
                        {
                            max_code=huffcodes[i];
                        }
                    }
                    i=0;
                    while(max_code&1)
                    {
                        i++;
                        max_code>>=1;
                    }
                    if(i>max_huff_len)
                    {
                        for(i=0; i<symbol_size; i++)
                        {
                            if(huffcodes[i]>max_code)
                            {
                                max_code=huffcodes[i];
                            }
                        }
                        printf("Error, %X>max huffcode(%i)\n", (int)max_code, max_huff_len);
                    }
                    if(max_code!=0)
                    {
                        for(i=0; i<symbol_size; i++)
                        {
                            if(huffcodes[i]>max_code)
                            {
                                max_code=huffcodes[i];
                            }
                        }
                        printf("Error rare maxcode: %X huffcode(%i)\n", (int)max_code, max_huff_len);
                    }
                }
            }
            pos+=delta;
            if(0)
            {
                symbol_count_t i;
                for(i=0; i<symbol_size; i++)
                {
                    if(s_len[i]!=0)
                    {
                        printf("%02lX, len=%i, huff=%lX\n", i, s_len[i], huffcodes[i]);
                    }
                }
            }

        }
    } while(delta>1);
    return totaal;
}

int main(int argc, char* argv[])
{
	uint32_t crc_table[256];
	int i;
	int max_huff_len;
	uint64_t start=rdtsc();
	uint64_t delta=4096;
	for(max_huff_len=16; max_huff_len<=16; max_huff_len++)
	{
	    i=1;
	    make_crc32_table(crc_table);
		uint64_t totaal=0;
		while(i<argc)
	    {
			uint32_t crc;
			int64_t size;
			uint8_t* data;
			const int symbol_size=256;
//			printf("%s: ",argv[i]);
			size=load_file(argv[i], &data);
			if(data==NULL)
		    {
				printf("Error loading %s\n", argv[i]);
				return -1;
			}
		    {
		        uint64_t result;
				result=encode(data, size, symbol_size, max_huff_len, delta);
//				printf("hufflen=%02i, size=%lu\n", max_huff_len, result);
				totaal+=result;
			}
			free(data);
			i++;
		}
		printf("Totaal(%02i) = %lu\n", max_huff_len, totaal);
	}
	printf("Tijd(delta=%li) totaal=%li make_huffmantable=%li, sort=%li, walktree=%li\n", delta, rdtsc()-start, time_huffman, time_sort, time_walktree);
    return 0;
}
