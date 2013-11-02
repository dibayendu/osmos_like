/* Texture loader for Linux using gdk-pixbuf-2.0 (a component of GTK).
 * Can load PNG, JPEG, GIF, BMP, TIFF, ..., etc.  
 *
 * Make sure the image has dimensions that are powers of two! (e.g., 32x64 is
 * ok, but 30x70 is not).
 *
 * Compile this file on Linux with 
 *  
 *   gcc -c texture_gdk.c `pkg-config gdk-pixbuf-2.0 --cflags --libs`
 *
 * or use the Makefile from the I3D website.
 */

#include <gdk-pixbuf/gdk-pixbuf.h>

#include "texture.h"

int is_power_2(int val)
{
    int count = 0;
    while (val)
    {
        count += val & 1;
        val >>= 1;
    }
    return 1;//count == 1;
}

int texture_is_valid_dimensions(int width, int height)
{
    return is_power_2(width) && is_power_2(height);
}

void flip_data(char *data, int pitch, int height)
{
    /* Flip the rows of the image data in-place */

    char *row1 = data;
    char *row2 = data + (height - 1) * pitch;
    int x, y;
    char tmp;

    for (y = 0; y < height >> 1; y++)
    {
        for (x = 0; x < pitch; x++)
        {
            tmp = row1[x];
            row1[x] = row2[x];
            row2[x] = tmp;
        }
        row1 += pitch;
        row2 -= pitch;
    }
}

GLuint texture_load_data(unsigned char *data, int width, int height, 
                         int components, int pitch,
                         GLint internalFormat, GLenum format, GLenum type)
{
    GLuint id;
    int alignment;
    int row_length;

    /* If pitch is negative, flip order of rows from top-to-bottom to
       bottom-to-top. */
    if (pitch < 0)
    {
        pitch = -pitch;
        flip_data(data, pitch, height);
    }

    if (pitch & 0x1)
        alignment = 1;
    else if (pitch & 0x2)
        alignment = 2;
    else
        alignment = 4;
    row_length = pitch / components;

    glPushClientAttrib(GL_CLIENT_PIXEL_STORE_BIT);
    glPixelStorei(GL_UNPACK_ALIGNMENT, alignment);
    glPixelStorei(GL_UNPACK_ROW_LENGTH, row_length);

    glGenTextures(1, &id);
    glBindTexture(GL_TEXTURE_2D, id);
    glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, width, height, 0, format, 
                 type, data);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFinish();

    glPopClientAttrib();

    return id;
}

GLuint texture_load(const char *filename)
{
    GdkPixbuf *pixbuf;
    int width, height, rowstride, channels;
    guchar *pixels;

    GLuint id;
    GLint internalformat;
    GLenum format;

    g_type_init();
    pixbuf = gdk_pixbuf_new_from_file(filename, NULL);

    if (!pixbuf)
    {
        fprintf(stderr, "Could not load file %s\n", filename);
        return 0;
    }

    width = gdk_pixbuf_get_width(pixbuf);
    height = gdk_pixbuf_get_height(pixbuf);

    if (!texture_is_valid_dimensions(width, height))
    {
        fprintf(stderr, "Image %s does not have dimensions of power 2\n", 
                filename);
        gdk_pixbuf_unref(pixbuf);
        return 0;
    }

    rowstride = gdk_pixbuf_get_rowstride(pixbuf);
    channels = gdk_pixbuf_get_n_channels(pixbuf);
    pixels = gdk_pixbuf_get_pixels(pixbuf);

    switch (channels) 
    {
        case 1:
            internalformat = format = GL_LUMINANCE;
            break;
        case 2:
            internalformat = format = GL_LUMINANCE_ALPHA;
            break;
        case 3:
            internalformat = format = GL_RGB;
            break;
        case 4:
            internalformat = format = GL_RGBA;
            break;
    }

    id = texture_load_data(pixels, width, height, 
                           channels, -rowstride, internalformat, 
                           format, GL_UNSIGNED_BYTE);

    gdk_pixbuf_unref(pixbuf);

    return id;
}
