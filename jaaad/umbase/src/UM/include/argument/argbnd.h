! ARGBND Control data calculated from NAMELIST-
     &  NBOUND_LOOKUP,                                                  &
      ! Headers from atmosphere boundary data sets
#if defined(ATMOS) && !defined(GLOBAL)
     &  FIXHD_BOUNDA,INTHD_BOUNDA,LOOKUP_BOUNDA,LOOKUP_COMP_BOUNDA,     &
     &  REALHD_BOUNDA,                                                  &
#endif
      ! Headers from ocean boundary data sets
#if (defined(OCEAN)&&defined(BOUNDSO))||(defined(OCEAN)&&defined(FLOOR))
     &  FIXHD_BOUNDO,INTHD_BOUNDO,LOOKUP_BOUNDO,REALHD_BOUNDO,          &
     &  joc_bdy_tracer,o_bdy_item_codes,o_bdy_prev,o_bdy_next,          &
#endif
