typedef struct
{
    FILE *fd;
    char* buffer;
    unsigned int len;
    unsigned int pos;
    unsigned int namelen;
} iohandle_s;
typedef iohandle_s* iohandle;

int SAVE_MD5_HASH = FALSE;
unsigned int SAVE_MD5_SIZE = 0;
int SAVE_MD5_SUCCESS = FALSE;

int CHECK_MD5_HASH = FALSE;

#define io_error(fd) fd == NULL

iohandle io_open(const char* filename)
{
    iohandle handle;

    if (CHECK_MD5_HASH)
    {
        FILE* fd = fopen(filename, "rb");
        if (!fd)
            return NULL;
        fseek(fd, 0L, SEEK_END);
        int filelen = ftell(fd);
        handle = malloc(sizeof(iohandle_s));
        handle->fd = NULL;
        handle->namelen = strlen(filename);
        handle->pos = handle->namelen + 5;
        handle->len = filelen + handle->pos;
        handle->buffer = malloc(handle->len);
        memcpy(handle->buffer, filename, handle->namelen);
        strcpy(handle->buffer + handle->namelen, ".md5");

        fseek(fd, 0L, SEEK_SET);
        filelen = fread(handle->buffer + handle->pos, 1, handle->len - handle->pos, fd);
        fclose(fd);
        fd = fopen(handle->buffer, "rb");
        if (!fd || filelen != handle->len - handle->pos)
        {
            remove(filename);
            free(handle->buffer);
            free(handle);
            return NULL;
        }

        char md5file[33];
        char md5data[33];
        fread(md5file, 1, 32, fd);
        fclose(fd);
        md5_raw_input(md5data, handle->buffer + handle->pos, handle->len - handle->pos);
        if (strncmp(md5file, md5data, 32))
        {
            remove(filename);
            remove(handle->buffer);
            free(handle->buffer);
            free(handle);
            return NULL;
        }

        return handle;
    }
    else
    {
        FILE* fd = fopen(filename, "rb");
        if (!fd)
            return NULL;
        handle = malloc(sizeof(iohandle_s));
        handle->fd = fd;
        handle->buffer = NULL;
        handle->len = 0;
        handle->pos = 0;
        return handle;
    }
}

iohandle io_create(const char *filename)
{
    iohandle handle;

    if (SAVE_MD5_HASH)
    {
        handle = malloc(sizeof(iohandle_s));
        handle->fd = NULL;
        handle->namelen = strlen(filename);
        handle->pos = handle->namelen + 5;
        handle->len = SAVE_MD5_SIZE + handle->pos;
        handle->buffer = malloc(handle->len);
        memcpy(handle->buffer, filename, handle->namelen);
        strcpy(handle->buffer + handle->namelen, ".md5");
        SAVE_MD5_SUCCESS = TRUE;
        return handle;
    }
    else
    {
        FILE* fd = fopen(filename, "wb");
        if (!fd)
            return NULL;
        handle = malloc(sizeof(iohandle_s));
        handle->fd = fd;
        handle->buffer = NULL;
        handle->len = 0;
        handle->pos = 0;
        return handle;
    }
}

int io_read(iohandle handle, void *buffer, unsigned int count)
{
    if (handle->fd != NULL)
        return fread(buffer, 1, count, handle->fd);
    if (handle->buffer == NULL)
        return 0;
    if (count > handle->len - handle->pos)
        count = handle->len - handle->pos;
    memcpy(buffer, handle->buffer + handle->pos, count);
    handle->pos += count;
    return count;
}

int io_write(iohandle handle, const void *buffer, unsigned int count)
{
    if (handle->fd != NULL)
        return fwrite(buffer, 1, count, handle->fd);
    if (handle->buffer == NULL || handle->len < handle->pos + count)
    {
        unsigned int newlen = handle->len*2;
        if (newlen < handle->pos + count)
            newlen = handle->pos + count;
        char *newbuf = malloc(newlen);
        if (handle->buffer != NULL && handle->pos > 0)
            memcpy(newbuf, handle->buffer, handle->pos);
        if (handle->buffer != NULL)
            free(handle->buffer);
        handle->buffer = newbuf;
        handle->len = newlen;
    }
    memcpy(handle->buffer + handle->pos, buffer, count);
    handle->pos += count;
    return count;
}

void writeThrough(char *filename, const void *buffer, unsigned int count)
{
#ifdef dwedwe_WIN32
    HANDLE hFile = CreateFileA(filename, GENERIC_WRITE, FILE_SHARE_READ, NULL, CREATE_ALWAYS, FILE_FLAG_WRITE_THROUGH, NULL);
    if (hFile == INVALID_HANDLE_VALUE)
    {
        SAVE_MD5_SUCCESS = FALSE;
        return;
    }
    DWORD written;
    if (!WriteFile(hFile, buffer, count, &written, NULL) || written != count)
        SAVE_MD5_SUCCESS = FALSE;
    CloseHandle(hFile);
#else
    FILE *fd = fopen(filename, "wb");
    if (!fd)
    {
        SAVE_MD5_SUCCESS = FALSE;
        return;
    }
    if (fwrite(buffer, 1, count, fd) != count)
        SAVE_MD5_SUCCESS = FALSE;
    fclose(fd);
#endif
}

void io_commit(iohandle handle)
{
    if (handle->fd != NULL)
    {
        fflush(handle->fd);
        return;
    }
    if (handle->buffer != NULL)
    {
        char md5hash[33];
        md5_raw_input(md5hash, handle->buffer + handle->namelen + 5, handle->pos - handle->namelen - 5);
        writeThrough(handle->buffer, md5hash, 32);
        handle->buffer[handle->namelen] = 0;
        writeThrough(handle->buffer, handle->buffer + handle->namelen + 5, handle->pos - handle->namelen - 5);
    }
}

void io_close(iohandle handle)
{
    if (handle->fd != NULL)
        fclose(handle->fd);
    if (handle->buffer != NULL)
        free(handle->buffer);
    free(handle);
}

#define io_reset remove
#define io_report(a,b,c,d,e,f,g)


/*void tempFileName (
    char	*buf, char c, giant NN)
{
    int remainder;

    remainder = gmodi(19999981, NN);
    sprintf (buf, "%c%07i", c, remainder % 10000000);
} */

int fileExists(char *filename)
{
    iohandle fd;
    fd = io_open(filename);
    if (io_error(fd))
        return FALSE;
    io_close(fd);
    return TRUE;
}

int io_read_uint32(iohandle fd, uint32_t* val, uint32_t* sum)
{
    if (io_read(fd, val, sizeof(uint32_t)) != sizeof(uint32_t))
        return FALSE;
    *sum += *val;
    return TRUE;
}

int io_write_uint32(iohandle fd, uint32_t val, uint32_t* sum)
{
    if (io_write(fd, &val, sizeof(uint32_t)) != sizeof(uint32_t))
        return FALSE;
    *sum += val;
    return TRUE;
}

int io_read_giant(iohandle fd, giant g, uint32_t *sum)
{
    uint32_t i, len, bytes;
    giant tmp;

    if (!io_read_uint32(fd, &len, sum))
        return FALSE;
    if (len == 0)
    {
        if (g != NULL)
            g->sign = 0;
        return TRUE;
    }
    if (gwdata != NULL)
    {
        if (len > ((uint32_t)gwdata->bit_length >> 5) + 5)
            return FALSE;
        tmp = g != NULL ? g : getg();
    }
    else
    {
        N = allocgiant(len + 10);
        tmp = N;
    }
    bytes = len*sizeof(uint32_t);
    if (io_read(fd, tmp->n, bytes) != bytes)
        return FALSE;
    tmp->sign = len;
    for (i = 0; i < len; i++)
        *sum += tmp->n[i];
    if (g == NULL && gwdata != NULL)
        freeg();
    return TRUE;
}

int io_write_giant(iohandle fd, giant g, uint32_t *sum)
{
    uint32_t i, len, bytes;

    len = g->sign;
    if (!io_write_uint32(fd, len, sum))
        return FALSE;
    bytes = len*sizeof(uint32_t);
    if (io_write(fd, g->n, bytes) != bytes)
        return FALSE;
    for (i = 0; i < len; i++)
        *sum += g->n[i];
    return TRUE;
}

int writeToFile(char *filename, uint32_t fingerprint, int32_t j, giant x, giant y)
{
    iohandle fd;
    uint32_t magicnum, version;
    uint32_t sum = 0;

    fd = io_create(filename);
    if (io_error(fd))
        return FALSE;

    magicnum = 0x9f2b3cd4;
    if (io_write(fd, &magicnum, sizeof(uint32_t)) != sizeof(uint32_t))
        goto writeerr;
    version = 0x00010001;
    if (io_write(fd, &version, sizeof(uint32_t)) != sizeof(uint32_t))
        goto writeerr;

    if (!io_write_uint32(fd, fingerprint, &sum))
        goto writeerr;
    if (!io_write_uint32(fd, (uint32_t)j, &sum))
        goto writeerr;

    if (!io_write_giant(fd, x, &sum))
        goto writeerr;
    if (y != NULL && !io_write_giant(fd, y, &sum))
        goto writeerr;
    if (y == NULL && !io_write_uint32(fd, 0, &sum))
        goto writeerr;

    if (io_write(fd, &sum, sizeof(uint32_t)) != sizeof(uint32_t))
        goto writeerr;

    io_commit(fd);
    io_close(fd);

    return TRUE;

writeerr:
    io_close(fd);
    return FALSE;
}

int writeToFileSafe(char* filename, uint32_t fingerprint, int32_t j, giant x, giant y)
{
    char newfilename[64];
    strcpy(newfilename, filename);
    strcat(newfilename, ".new");

    if (!writeToFile(newfilename, fingerprint, j, x, y))
    {
        remove(newfilename);
        return FALSE;
    }

    remove(filename);
    rename(newfilename, filename);

    return TRUE;
}

int writeToFileMD5(char *filename, uint32_t fingerprint, int32_t j, giant x, giant y)
{
    int retval = TRUE;

    SAVE_MD5_HASH = TRUE;
    SAVE_MD5_SIZE = ((int)gwdata->bit_length >> 3)*(y != NULL ? 2 : 1) + 128;
    if ((retval = writeToFile(filename, fingerprint, j, x, y)))
        retval = SAVE_MD5_SUCCESS;
    SAVE_MD5_HASH = FALSE;

    return retval;
}

int readFromFile(char *filename, uint32_t fingerprint, int32_t *j, giant x, giant y)
{
    iohandle fd;
    uint32_t magicnum, version;
    uint32_t i, sum = 0;

    fd = io_open(filename);
    if (io_error(fd))
        return FALSE;

    if (io_read(fd, &magicnum, sizeof(uint32_t)) != sizeof(uint32_t))
        goto readerr;
    if (magicnum != 0x9f2b3cd4)
        goto readerr;

    if (io_read(fd, &version, sizeof(uint32_t)) != sizeof(uint32_t))
        goto readerr;
    if (version != 0x00010001)
        goto readerr;

    if (!io_read_uint32(fd, &i, &sum))
        goto readerr;
    if (i != fingerprint)
        goto readerr;

    *j = 0;
    if (!io_read_uint32(fd, (uint32_t*)j, &sum))
        goto readerr;

    if (!io_read_giant(fd, x, &sum))
        goto readerr;
    if (!io_read_giant(fd, y, &sum))
        goto readerr;

    if (io_read(fd, &i, sizeof(uint32_t)) != sizeof(uint32_t))
        goto readerr;
    if (i != sum)
        goto readerr;

    io_close(fd);
    return TRUE;

readerr:
    io_close(fd);
    io_reset(filename);
    return FALSE;
}

int readFromFileMD5(char *filename, uint32_t fingerprint, int32_t* j, giant x, giant y)
{
    CHECK_MD5_HASH = TRUE;
    int retval = readFromFile(filename, fingerprint, j, x, y);
    CHECK_MD5_HASH = FALSE;
    return retval;
}

// Writes factor to factors.txt
void report_factor(giant f)
{
    char *buf = malloc(abs(f->sign)*10 + 10);
    gtoc(f, buf, abs(f->sign)*10 + 10);
    printf("Found factor %s\n", buf);
    FILE *fp = fopen("factors.txt", "a");
    if (fp)
    {
        fprintf(fp, "%s | %s\n", buf, Nstr);
        fclose(fp);
    }
}
