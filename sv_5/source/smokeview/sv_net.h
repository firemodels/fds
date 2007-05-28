#ifdef pp_SVNET
typedef struct _svcom {
  char label[256];
  unsigned int src_id, dest_id;
  unsigned char dest_ip[4];
  int command_type;
  char char_data[1024];
  float  rvals[16];
  int ivals[16];
  struct _svcom *prev, *next;
} svcom;

typedef struct _tcpdata {
  char label[256];
  char command[1024];
  unsigned int src_id, dest_id;
  unsigned char dest_ip[4];
  struct _tcpdata *prev, *next;
} tcpdata;

#endif