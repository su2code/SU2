c
c     file open modes
c
      integer CGIO_MODE_READ, CGIO_MODE_WRITE, CGIO_MODE_MODIFY
      parameter (CGIO_MODE_READ   = 0)
      parameter (CGIO_MODE_WRITE  = 1)
      parameter (CGIO_MODE_MODIFY = 2)
c
c     database file types
c
      integer CGIO_FILE_NONE, CGIO_FILE_ADF, CGIO_FILE_HDF5,
     &        CGIO_FILE_ADF2
      parameter (CGIO_FILE_NONE = 0)
      parameter (CGIO_FILE_ADF  = 1)
      parameter (CGIO_FILE_HDF5 = 2)
      parameter (CGIO_FILE_ADF2 = 3)
c
c     dimension limits
c
      integer CGIO_MAX_DATATYPE_LENGTH, CGIO_MAX_DIMENSIONS,
     &        CGIO_MAX_NAME_LENGTH, CGIO_MAX_LABEL_LENGTH,
     &        CGIO_MAX_VERSION_LENGTH, CGIO_MAX_ERROR_LENGTH,
     &        CGIO_MAX_LINK_DEPTH, CGIO_MAX_FILE_LENGTH,
     &        CGIO_MAX_LINK_LENGTH
      parameter (CGIO_MAX_DATATYPE_LENGTH = 2)
      parameter (CGIO_MAX_DIMENSIONS      = 12)
      parameter (CGIO_MAX_NAME_LENGTH     = 32)
      parameter (CGIO_MAX_LABEL_LENGTH    = 32)
      parameter (CGIO_MAX_VERSION_LENGTH  = 32)
      parameter (CGIO_MAX_ERROR_LENGTH    = 80)
      parameter (CGIO_MAX_LINK_DEPTH      = 100)
      parameter (CGIO_MAX_FILE_LENGTH     = 1024)
      parameter (CGIO_MAX_LINK_LENGTH     = 4096)

