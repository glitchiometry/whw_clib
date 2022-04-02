#include "system.h"

void mkdir_s(char *dirname)
{
	mkdir(dirname, S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH); // check these flags
}
