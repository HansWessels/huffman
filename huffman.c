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


#define MAX_SYMBOL_SIZE 256
#define MAX_HUFFMAN_LEN 16

#if MAX_SYMBOL_SIZE <= 256 /* voor max_symbol_size <=256 */
    typedef uint8_t symbol_t;
    #define PRI_SYMBOL_T PRIu8
    typedef int16_t symbol_count_t;
    #define PRI_SYMBOL_COUNT_T PRIi16
#elif MAX_SYMBOL_SIZE <= 32768 /* voor max_simbol_size <=32768 */
    typedef uint16_t symbol_t;
    #define PRI_SYMBOL_T PRIu16
    typedef int16_t symbol_count_t;
    #define PRI_SYMBOL_COUNT_T PRIi16
#elif MAX_SYMBOL_SIZE <= 65536 /* voor max_simbol_size <=65536 */
    typedef uint16_t symbol_t;
    #define PRI_SYMBOL_T PRIu16
    typedef int32_t symbol_count_t;
    #define PRI_SYMBOL_COUNT_T PRIi32
#elif MAX_SYMBOL_SIZE <= (1<<31) /* voor max_simbol_size <=1<<31 */
    typedef uint32_t symbol_t;
    #define FORMAT_SYMBOL_T PRIu32
    typedef int32_t symbol_count_t;
    #define PRI_SYMBOL_COUNT_T PRIi32
#elif MAX_SYMBOL_SIZE <= (1<<32) /* voor max_simbol_size <=1<<32 */
    typedef uint32_t symbol_t;
    #define PRI_SYMBOL_T PRIu32
    typedef int64_t symbol_count_t;
    #define F_SYMBOL_COUNT_T PRIi64
#elif 0 /* rest */
    typedef uint64_t symbol_t;
    #define PRI_SYMBOL_T PRIu64
    typedef int64_t symbol_count_t;
    #define PRI_SYMBOL_COUNT_T PRIi64
#endif


#if MAX_HUFFMAN_LEN <= 8
    typedef uint8_t huffman_t;
    #define PRI_HUFFMAN_T PRIX8
#elif MAX_HUFFMAN_LEN <= 16
    typedef uint16_t huffman_t;
    #define PRI_HUFFMAN_T PRIX16
#elif MAX_HUFFMAN_LEN <= 32
    typedef uintt32_t huffman_t;
    #define PRI_HUFFMAN_T PRIX32
#else
    typedef uint64_t huffman_t;
    #define PRI_HUFFMAN_T PRIX64
#endif
typedef uint64_t freq_t;

uint64_t time_huffman=0;
uint64_t time_sort=0;
uint64_t time_walktree=0;
uint64_t time_normal_tree=0;

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

void freq_count(const uint8_t *data, uint64_t size, freq_t freq[], symbol_count_t symbol_size)
{
    symbol_count_t i;
    uint64_t pos;
    for(i=0; i<symbol_size; i++)
    {
        freq[i]=0;
    }
    huffman_t symbol_mask=0;
    uint64_t current_value=0;
    symbol_size>>=1;
    while(symbol_size!=0)
    {
        symbol_size>>=1;
        symbol_mask<<=1;
        symbol_mask^=1;
    }
    current_value=(uint64_t)data[0]+((uint64_t)data[1]<<8)+((uint64_t)data[2]<<16)+((uint64_t)data[3]<<24)+((uint64_t)data[4]<<32)+((uint64_t)data[5]<<40)+((uint64_t)data[6]<<48);
    current_value<<=8;
    for(pos=0; pos<size;pos++)
    {
        symbol_t tmp;
        current_value>>=8;
        current_value+=(uint64_t)data[pos+7]<<56;
        tmp=(symbol_t)(current_value&symbol_mask);
        freq[tmp]++;
    }
}

typedef struct
{
    freq_t freq;
    symbol_t symbol;
} freq_compare_t;

int freq_compare(const void* a_in, const void* b_in)
{
    const freq_compare_t* a=(freq_compare_t *)a_in;
    const freq_compare_t* b=(freq_compare_t *)b_in;
    if (a->freq < b->freq)
    {
        return -1;
    }
    if (a->freq > b->freq)
    {
        return 1;
    }
    return 0;
}

/*
 * sorten met Q sort is ongeveer twee maal langzamer als radix sort
 */
void sort_symbols(symbol_t symbols[MAX_SYMBOL_SIZE], freq_t freq[MAX_SYMBOL_SIZE], symbol_count_t symbol_count)
{
    freq_compare_t pairs[MAX_SYMBOL_SIZE];
    symbol_count_t i;
    for(i=0; i<symbol_count; i++)
    {
        pairs[i].freq=freq[i];
        pairs[i].symbol=symbols[i];
    }
    qsort(pairs, symbol_count, sizeof(freq_compare_t), freq_compare);
    for(i=0; i<symbol_count; i++)
    {
        freq[i]=pairs[i].freq;
        symbols[i]=pairs[i].symbol;
    }
}

/*
 * Het lijkt er op dat het geen zin heeft om voor kleine datasets insertion sort te
 * gebruiken i.p.v. radix sort, ik heb geen significante tijd wints kunnen meten.
 */

 #define INSERTION_GRENS 16

void insertion_sort_symbols(symbol_t symbols[MAX_SYMBOL_SIZE], freq_t freq[MAX_SYMBOL_SIZE], symbol_count_t symbol_count)
{
    symbol_count_t i;
    for(i=1; i<symbol_count; i++)
    {
        freq_t cur_freq=freq[i];
        symbol_t cur_symbol=symbols[i];
        symbol_count_t pos=i-1;
        while(pos>=0 && (freq[pos]>cur_freq))
        {
            freq[pos+1]=freq[pos];
            symbols[pos+1]=symbols[pos];
            pos--;
        }
        freq[pos+1]=cur_freq;
        symbols[pos+1]=cur_symbol;
    }
}

void radix_sort_symbols(symbol_t symbols[MAX_SYMBOL_SIZE], freq_t freq[MAX_SYMBOL_SIZE], symbol_count_t symbol_count)
{
    symbol_t tmp_symbols[MAX_SYMBOL_SIZE];
    freq_t tmp_freq[MAX_SYMBOL_SIZE];

    #define BUCKET_BITS 8
    #define BUCKET_SIZE (1<<BUCKET_BITS)
    #define BUCKET_MASK (BUCKET_SIZE-1)

    symbol_count_t bucket[BUCKET_SIZE];
    freq_t max=0;
    int shift=0;
    if(symbol_count<INSERTION_GRENS)
    {
        insertion_sort_symbols(symbols, freq, symbol_count);
    }
    symbol_count_t i=symbol_count;
    do
    {
        i--;
        if(max<freq[i])
        {
            max=freq[i];
        }
    } while(i>0);

    for(;;)
    {
        symbol_count_t i;
        i=symbol_count;
        memset(bucket, 0, sizeof(bucket));
        do
        {
            i--;
            bucket[((freq[i]>>shift) & BUCKET_MASK)]++;
        } while(i>0);

        symbol_count_t start;
        start=0;
        for(i=0; i<BUCKET_SIZE; i++)
        {
            symbol_count_t tmp=bucket[i];
            bucket[i]=start;
            start+=tmp;
        }
        for(i=0; i<symbol_count; i++)
        {
            int tmp=(freq[i]>>shift)&BUCKET_MASK;
            tmp_freq[bucket[tmp]]=freq[i];
            tmp_symbols[bucket[tmp]]=symbols[i];
            bucket[tmp]++;
        }
        max>>=BUCKET_BITS;
        if(max==0)
        {
            memcpy(freq, tmp_freq, symbol_count*sizeof(freq[0]));
            memcpy(symbols, tmp_symbols, symbol_count*sizeof(symbols[0]));
            return;
        }
        shift+=BUCKET_BITS;
        i=symbol_count;
        memset(bucket, 0, sizeof(bucket));
        do
        {
            i--;
            bucket[((tmp_freq[i]>>shift) & BUCKET_MASK)]++;
        } while(i>0);
        start=0;
        for(i=0; i<BUCKET_SIZE; i++)
        {
            symbol_count_t tmp=bucket[i];
            bucket[i]=start;
            start+=tmp;
        }
        for(i=0; i<symbol_count; i++)
        {
            int tmp=(tmp_freq[i]>>shift)&BUCKET_MASK;
            freq[bucket[tmp]]=tmp_freq[i];
            symbols[bucket[tmp]]=tmp_symbols[i];
            bucket[tmp]++;
        }
        max>>=BUCKET_BITS;
        if(max==0)
        {
            return;
        }
        shift+=BUCKET_BITS;
    }
}

int huffman_sanety_check(int s_len[], huffman_t huff_codes[], symbol_count_t symbol_count, int max_huff_len)
{
    int max_len=0;
    symbol_count_t max_len_symbol;
    symbol_count_t i;
    int result=0;
    i=symbol_count;
    max_len_symbol=0;
    do
    {
        i--;
        if(s_len[i]>max_len)
        {
            max_len=s_len[i];
            max_len_symbol++;
        }
    } while(i>0);
    if((max_len_symbol<2) && (max_len<=1))
    { /* aantal symbolen <2, noets te checken */
        return result;
    }
    /* max_len > max_huff_len? */
    if(max_len>max_huff_len)
    {
        i=symbol_count;
        do
        {
            i--;
            if(s_len[i]>max_huff_len)
            {
                fprintf(stderr, "Te lange huffmancode: s_len[%"PRI_SYMBOL_COUNT_T"]=%i, huffcode=%" PRI_HUFFMAN_T " > %i (MAX_HUFF_LEN)\n", i, s_len[i], huff_codes[i], max_huff_len);
                result++;
            }
        }
        while(i>0);
    }
    /* is geen van de huffman_codes > max_len? */
    {
        uint64_t mask;
        mask=0;
        mask=~mask;
        mask<<=max_len;
        i=symbol_count;
        do
        {
            i--;
            if((huff_codes[i]&mask)!=0)
            {
                fprintf(stderr, "Te lange huffmancode: s_len[%"PRI_SYMBOL_COUNT_T"]=%i, huffcode=%" PRI_HUFFMAN_T "\n", i, s_len[i], huff_codes[i]);
                result++;
            }
        } while(i>0);
    }
    /* klopt de tree? */
    {
        int len_count[MAX_HUFFMAN_LEN+1]={0};
        huffman_t huffcode[MAX_HUFFMAN_LEN+1];
        huffman_t max_code;
        symbol_count_t i;
        i=symbol_count;
        do
        {
            i--;
            len_count[s_len[i]]++;
        } while(i>0);
        huffman_t start_huffcode=0;
        huffcode[0]=0;
        for(i=1; i<=MAX_HUFFMAN_LEN; i++)
        {
            start_huffcode<<=1;
            huffcode[i]=start_huffcode;
            start_huffcode+=len_count[i];
        }
        for(i=0; i<symbol_count; i++)
        {
            if(huff_codes[i]!=huffcode[s_len[i]])
            {
                fprintf(stderr, "Foute huffman code: huff_code[%"PRI_SYMBOL_COUNT_T"]=%" PRI_HUFFMAN_T "!=%" PRI_HUFFMAN_T "\n", i, huff_codes[i], huffcode[s_len[i]]);
                result++;
            }
            if(s_len[i]!=0)
            {
                huffcode[s_len[i]]++;
            }
        }
        i=symbol_count;
        max_code=0;
        do
        {
            i--;
            if(huff_codes[i]>max_code)
            {
                max_code=huff_codes[i];
                max_len_symbol=i;
            }
        } while(i>0);
        i=max_len;
        huffman_t tmp=max_code;
        do
        {
            i--;
            if((tmp&1)!=1)
            {
                fprintf(stderr, "Foute huffman maximale code, niet allemaal 1: huff_code[%"PRI_SYMBOL_COUNT_T"]=%" PRI_HUFFMAN_T "!=%" PRI_HUFFMAN_T "\n", max_len_symbol, huff_codes[max_len_symbol], ~((~(huffman_t)0)<<max_len));
                result++;
            }
            tmp>>=1;
        } while(i>0);
        if(tmp!=0)
        {
            fprintf(stderr, "Foute huffman maximale code, te veel 1: huff_code[%"PRI_SYMBOL_COUNT_T"]=%" PRI_HUFFMAN_T "!=%" PRI_HUFFMAN_T "\n", max_len_symbol, huff_codes[max_len_symbol], ~((~(huffman_t)0)<<max_len));
            result++;
        }
    }
    return result;
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
    huffcode[0]=0;
    for(i=1; i<=MAX_HUFFMAN_LEN; i++)
    {
        start_huffcode<<=1;
        huffcode[i]=start_huffcode;
        start_huffcode+=len_count[i];
    }
    for(i=0; i<symbol_count; i++)
    {
        huff_codes[i]=huffcode[s_len[i]];
        huffcode[s_len[i]]+=(s_len[i]!=0);
    }
}

int make_huffman_table(int s_len[], huffman_t huff_codes[], const freq_t in_freq[], int max_huff_len, const symbol_count_t symbol_size)
{
    static freq_t freq_array[MAX_SYMBOL_SIZE*2];
    static symbol_count_t tree_array[MAX_HUFFMAN_LEN*MAX_SYMBOL_SIZE*2];
    static symbol_t symbols[MAX_SYMBOL_SIZE];
    static int len[2*MAX_SYMBOL_SIZE];
    symbol_count_t* tree=tree_array;
    freq_t* freq=freq_array;
//    freq_t* pairs_freq=freq+MAX_SYMBOL_SIZE;
    symbol_count_t symbol_count;
    symbol_count_t pairs_count;
    symbol_count_t i;
    const uint64_t een=1;
    if(max_huff_len>MAX_HUFFMAN_LEN)
    {
        fprintf(stderr, "max_huff_len te groot: %i\n", max_huff_len);
        return -1;
    }
    if(max_huff_len<64)
    {
        if((een<<max_huff_len)<symbol_size)
        {
            fprintf(stderr, "(1<<max_huff_len) < symbol_size : 1<<%" PRIi64 " < %" PRI_SYMBOL_COUNT_T "\n", een<<max_huff_len, symbol_size);
            return -1;
        }
    }
    if(symbol_size>MAX_SYMBOL_SIZE)
    {
        fprintf(stderr, "symbolsize(%i)>MAX_SYMBOL_SIZE(%i)\n", symbol_size, MAX_SYMBOL_SIZE);
        return -1;
    }
    i=symbol_size;
    symbol_count=0;
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
    { /* special cases 0, 1 en 2 symbolen */
        huffman_t code=0;
        memset(s_len, 0, symbol_size*sizeof(s_len[0]));
        memset(huff_codes, 0, symbol_size*sizeof(huff_codes[0]));
        for(i=0; i<symbol_size; i++)
        {
            if(in_freq[i]!=0)
            {
                huff_codes[i]=code;
                code++;
                s_len[i]=1;
            }
        }
        return 0;
    }
    {
        uint64_t time_start=rdtsc();
        radix_sort_symbols(symbols, freq, symbol_count);
        //sort_symbols(symbols, freq, symbol_count);
        time_sort+=rdtsc()-time_start;
    }
    if(0)
    { /* eerst een traditionele huffmanboom bouwen */
        uint64_t time_start=rdtsc();
        freq_t parent_freq;
        symbol_count_t tree_pos=0;
        symbol_count_t symbol_pos=0;
        symbol_count_t pair_pos=1;
        symbol_count_t node;
        freq_t* pairs_freq=freq+MAX_SYMBOL_SIZE;
        freq[symbol_count]=~0; /* sentry */
        for(node=1; node<symbol_count; node++)
        {
            pairs_freq[node]=~0; /* sentries */
        }
        node=1;
        do
        {
            if(freq[symbol_pos]<=pairs_freq[pair_pos])
            {
                tree[tree_pos]=symbol_pos;
                parent_freq=freq[symbol_pos];
                symbol_pos++;
            }
            else
            {
                tree[tree_pos]=-pair_pos;
                parent_freq=pairs_freq[pair_pos];
                pair_pos++;
            }
            tree_pos++;
            if(freq[symbol_pos]<=pairs_freq[pair_pos])
            {
                tree[tree_pos]=symbol_pos;
                parent_freq+=freq[symbol_pos];
                symbol_pos++;
            }
            else
            {
                tree[tree_pos]=-pair_pos;
                parent_freq+=pairs_freq[pair_pos];
                pair_pos++;
            }
            tree_pos++;
            tree[tree_pos]=node;
            pairs_freq[node]=parent_freq;
            node++;
            tree_pos++;
        } while(tree_pos>(3*symbol_count));

        { /* bouw s_len */
            freq[-node]=0;
            do
            {
                tree_pos--;
                parent_freq=freq[tree[tree_pos]];
                tree_pos--;
                freq[tree[tree_pos]]=parent_freq+1;
                tree_pos--;
                freq[tree[tree_pos]]=parent_freq+1;
            } while(tree_pos>0);
            if(parent_freq<max_huff_len)
            {
                for(i=0; i<symbol_count; i++)
                {
                    s_len[symbols[i]]=(int)freq[i];
                }
                make_huffman_codes(s_len, huff_codes, symbol_size);
                time_normal_tree+=rdtsc()-time_start;
                return 0;
            }
        }
        time_normal_tree+=rdtsc()-time_start;
    }
    freq+=symbol_count; /* freq array index 1 based */
    freq[0]=0; /* sentry */
    {
        symbol_count_t symbol_pos;
        pairs_count=symbol_count>>1;
        symbol_pos=-1-(symbol_count&1);
        i=pairs_count;
        do
        { /* eerste ronde, gewoon de gegeven character bij elkaar voegen */
            i--;
            freq_t node_freq;
            node_freq=freq[symbol_pos];
            tree[i*2+1]=symbol_pos;
            symbol_pos--;
            node_freq+=freq[symbol_pos];
            tree[i*2]=symbol_pos;
            symbol_pos--;
            freq[i+1]=node_freq;
        } while(i>0);
        max_huff_len--;
    }
    do
    { /* merge symbols, max_huff_len-1 keer */
        symbol_count_t symbol_pos=-1;
        symbol_count_t pair_pos=pairs_count;
        freq_t next_symbol_freq=freq[symbol_pos];
        freq_t next_pair_freq=freq[pair_pos];
        tree+=symbol_count*2;
        pairs_count=symbol_count+pairs_count;
        if((pairs_count)&1)
        { /* oneven som, waarde met hoogste freq doet niet mee */
            if(next_pair_freq>next_symbol_freq)
            {
                pair_pos--;
                next_pair_freq=freq[pair_pos];
            }
            else
            {
                symbol_pos--;
                next_symbol_freq=freq[symbol_pos];
            }
        }
        pairs_count>>=1;
        i=pairs_count;
        do
        { /* maak de nieuwe pairs */
            freq_t node_freq;
            i--;
            if(next_pair_freq>next_symbol_freq)
            {
                node_freq=next_pair_freq;
                tree[i*2+1]=0;
                pair_pos--;
                next_pair_freq=freq[pair_pos];
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
                tree[i*2]=0;
                pair_pos--;
                next_pair_freq=freq[pair_pos];
            }
            else
            {
                node_freq+=next_symbol_freq;
                tree[i*2]=symbol_pos;
                symbol_pos--;
                next_symbol_freq=freq[symbol_pos];
            }
            freq[i+1]=node_freq;
        } while(i>0);
        max_huff_len--;
    } while(max_huff_len>0);
    memset(s_len, 0, symbol_size*sizeof(s_len[0]));
    {
        uint64_t time_start=rdtsc();
        symbol_count_t i;
        memset(freq-symbol_count, 0, symbol_count*(sizeof(freq[0])));
        i=2*(symbol_count-1);  /* N symbolen levert altijd N-1 pairs op */
        do
        {
            freq[0]=0;
            do
            {
                i--;
                freq[tree[i]]++;
            } while(i>0);
            i=2*freq[0];
            tree-=symbol_count*2;
        } while(i>0);
        i=symbol_count;
        freq-=symbol_count; /* undo negatieve symbol index */
        do
        {
            i--;
            s_len[symbols[i]]=(int)freq[i];
        } while(i>0);
        time_walktree+=rdtsc()-time_start;
    }
    make_huffman_codes(s_len, huff_codes, symbol_size);
    return 0;
}

uint64_t encode(const uint8_t *data, const int64_t size, const symbol_count_t symbol_size, const int max_huff_len, const uint64_t block_size)
{
    freq_t freq[MAX_SYMBOL_SIZE]={0};
    int s_len[MAX_SYMBOL_SIZE];
    huffman_t huffcodes[MAX_SYMBOL_SIZE];
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
                huffman_sanety_check(s_len, huffcodes, symbol_size, max_huff_len);
            }
            pos+=delta;
            if(0)
            {
                symbol_count_t i;
                for(i=0; i<symbol_size; i++)
                {
                    if(s_len[i]!=0)
                    {
                        printf("%" PRI_SYMBOL_COUNT_T", len=%i, huff=%"PRI_HUFFMAN_T"\n", i, s_len[i], huffcodes[i]);
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
	uint64_t start;
	uint64_t delta=4096000000;
	int symbol_size;
	symbol_size=256;
	do
	{
	int symbol_bits=0;
	int tmp=symbol_size>>1;
	while(tmp!=0)
	{
	    symbol_bits++;
		tmp>>=1;
	}
	max_huff_len=symbol_bits+8;
	for(max_huff_len=symbol_bits; max_huff_len<=symbol_bits+8; max_huff_len++)
	{
	    i=1;
	    make_crc32_table(crc_table);
		start=rdtsc();
		time_huffman=0;
		time_sort=0;
		time_walktree=0;
		uint64_t totaal=0;
		while(i<argc)
	    {
			uint32_t crc;
			int64_t size;
			uint8_t* data;
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
//        printf("Tijd(delta=%li, symbol_size=%i) totaal=%li make_huffmantable=%li, sort=%li, walktree=%li\n", delta, symbol_size, rdtsc()-start, time_huffman, time_sort, time_walktree);
	}
	symbol_size>>=1;
	} while(0); //while(symbol_size>1);
    return 0;
}
