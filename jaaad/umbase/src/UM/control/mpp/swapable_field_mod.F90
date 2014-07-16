    MODULE Swapable_Field_Mod
      TYPE swapable_field_pointer_type
        INTEGER :: field_type
        INTEGER :: levels
        INTEGER :: rows
        LOGICAL :: vector
        REAL, POINTER :: field(:, :, :)
      END TYPE swapable_field_pointer_type
    END MODULE Swapable_Field_Mod
