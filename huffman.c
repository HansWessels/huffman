/*
** testen van een Huffman encoder
** in het bijzonder het package merge algoritme
*/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SYMBOL_SIZE 512
#define SYMBOL_SIZE 256
#define MAX_HUFFMAN_LEN 64

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

void freq_count(const uint8_t *data, int64_t size, uint64_t freq[])
{
	while(size>0)
	{
		size--;
		freq[data[size]]++;
	}
}

int freq_compare(const void* a_in, const void* b_in)
{
    uint64_t* a=(uint64_t *)a_in;
    uint64_t* b=(uint64_t *)b_in;
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

void sort_symbols(int symbols[MAX_SYMBOL_SIZE], uint64_t freq[MAX_SYMBOL_SIZE], int symbol_count)
{
    uint64_t pairs[2*MAX_SYMBOL_SIZE];
    int i;
    for(i=0; i<symbol_count; i++)
    {
        pairs[2*i]=freq[symbols[i]];
        pairs[2*i+1]=symbols[i];
    }
    qsort(pairs, symbol_count, 2*sizeof(uint64_t), freq_compare);
    for(i=0; i<symbol_count; i++)
    {
        symbols[i]=(int)pairs[2*i+1];
    }
}

void walk_tree(int node, int c_left[MAX_SYMBOL_SIZE], int c_right[MAX_SYMBOL_SIZE], int s_len[])
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

void make_huffman_codes(int s_len[], uint64_t huff_codes[], int symbol_count)
{
    int len_count[MAX_HUFFMAN_LEN+1]={0};
    uint64_t huffcode[MAX_HUFFMAN_LEN+1];
    int i;
    i=symbol_count;
    do
    {
        i--;
        len_count[s_len[i]]++;
    } while(i>0);
    uint64_t start_huffcode=0;
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

int make_huffman_table(int s_len[], uint64_t huff_codes[], const uint64_t in_freq[], const int max_huff_len, const int symbol_size)
{
    uint64_t freq_array[MAX_SYMBOL_SIZE+MAX_SYMBOL_SIZE*MAX_HUFFMAN_LEN]={0};
    int c_left_array[MAX_SYMBOL_SIZE*MAX_HUFFMAN_LEN];
    int c_right_array[MAX_SYMBOL_SIZE*MAX_HUFFMAN_LEN];
    int symbols[MAX_SYMBOL_SIZE];
    int pairs[2*MAX_SYMBOL_SIZE];
    uint64_t* freq=freq_array+MAX_SYMBOL_SIZE*MAX_HUFFMAN_LEN;
    int* c_left=c_left_array+MAX_SYMBOL_SIZE*MAX_HUFFMAN_LEN;
    int* c_right=c_right_array+MAX_SYMBOL_SIZE*MAX_HUFFMAN_LEN;
    int symbol_count=0;
    int pairs_count;
    int node=-1;
    int i;
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
            printf("(1<<max_huff_len) < symbol_size : 1<<%i < %i", max_huff_len, symbol_size);
            return -1;
        }
    }
    i=symbol_size;
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
    i=max_huff_len;
    pairs_count=0;
    do
    { /* merge symbols, max_huff_len keer */
        int merge_pos=pairs_count+symbol_count-1;
        int symbol_pos=symbol_count-1;
        int pair_pos=pairs_count-1;
        do
        { /* merge symbols met pairs */
            if(pair_pos<0)
            {
                pairs[merge_pos]=symbols[symbol_pos];
                merge_pos--;
                symbol_pos--;
            }
            else if(freq[symbols[symbol_pos]]>=freq[pairs[pair_pos]])
            {
                pairs[merge_pos]=symbols[symbol_pos];
                merge_pos--;
                symbol_pos--;
            }
            else
            {
                pairs[merge_pos]=pairs[pair_pos];
                merge_pos--;
                pair_pos--;
            }
        } while(symbol_pos>=0);
        int total_pairs=(pairs_count+symbol_count)/2;
        pairs_count=0;
        do
        { /* maak de nieuwe pairs */
            freq[node]=freq[pairs[pairs_count*2]]+freq[pairs[pairs_count*2+1]];
            c_left[node]=pairs[pairs_count*2];
            c_right[node]=pairs[pairs_count*2+1];
            pairs[pairs_count]=node;
            pairs_count++;
            total_pairs--;
            node--;
        } while(total_pairs>0);
        i--;
    } while(i>0);
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

int encode(const uint8_t *data, const int64_t size, const int symbol_size, const int max_huff_len)
{
    uint64_t freq[SYMBOL_SIZE]={0};
    int s_len[SYMBOL_SIZE];
    uint64_t huffcodes[SYMBOL_SIZE];
    int i;
    freq_count(data, size, freq);
    make_huffman_table(s_len, huffcodes, freq, max_huff_len, symbol_size);
    {
        uint64_t totaal=0;
        int i;
        for(i=0; i<symbol_size; i++)
        {
            totaal+=freq[i]*s_len[i];
        }
        printf("Total_bits(max_huf=%i)=%lu\n", max_huff_len, totaal);
    }
    return 0;
}

int main(int argc, char* argv[])
{
	uint32_t crc_table[256];
	int i;
	i=1;
	make_crc32_table(crc_table);
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
		encode(data, size, symbol_size, max_huff_len);
		free(data);
		i++;
	}
    return 0;
}
