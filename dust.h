#ifndef DUST_H
#define DUST_H


#ifdef __cplusplus
extern "C" {
#endif

#define MAXREG    1001

  typedef struct {
    int from;
    int to;
    int score;
  } REGION;
  
  static REGION reg[MAXREG];
  static int nreg;

  static int word = 3; 
  static int window = 64; 
  static int window2 = 32; 
  static int level = 20;
  
  static int mv, iv, jv;
  
  void set_dust_level(int value);
  
  void set_dust_window(int value);
  
  void set_dust_word(int value);

  static void wo1(int len, char *s, int ivv);

  static int wo(int len, char *s, int *beg, int *end);

  void dust(int len, char *s, int *sum);
  
  REGION *dust_segs(int len, char *s);
  
  REGION *dust_mix(int len, char *s);
  
#ifdef __cplusplus
}
#endif

#endif
