  E!  n   k820309              10.1        �ǰM                                                                                                       
       lensing.f90 LENSING          LENS_CLS LENSING_INCLUDES_TENSORS LENSING_METHOD LENSING_METHOD_FLAT_CORR LENSING_METHOD_CURV_CORR LENSING_METHOD_HARMONIC BESSI BESSJ0                                        
                                              
                    �                         
                 @  @                        '�     0      #WANTCLS    #WANTTRANSFER    #WANTSCALARS    #WANTTENSORS    #WANTVECTORS 	   #DOLENSING 
   #WANT_ZSTAR    #WANT_ZDRAG    #NONLINEAR    #MAX_L    #MAX_L_TENSOR    #MAX_ETA_K    #MAX_ETA_K_TENSOR    #OMEGAB    #OMEGAC    #OMEGAV    #OMEGAN    #H0    #TCMB    #YHE    #NUM_NU_MASSLESS    #NUM_NU_MASSIVE    #NU_MASS_SPLITTINGS    #NU_MASS_EIGENSTATES    #NU_MASS_DEGENERACIES    #NU_MASS_FRACTIONS    #SCALAR_INITIAL_CONDITION    #OUTPUTNORMALIZATION     #ACCURATEPOLARIZATION !   #ACCURATEBB "   #ACCURATEREIONIZATION #   #MASSIVENUMETHOD $   #INITPOWER %   #REION /   #RECOMB 7   #TRANSFER =   #INITIALCONDITIONVECTOR D   #ONLYTRANSFERS E   #REIONHIST F   #FLAT N   #CLOSED O   #OPEN P   #OMEGAK Q   #CURV R   #R S   #KSIGN T   #TAU0 U   #CHI0 V           � D                                             � D                                            � D                                            � D                                            � D                      	                      � D                      
                      � D                                            � D                                            � D                               	              � D                           $   
              � D                           (                 � D                          0      
           � D                          8      
           � D                          @      
           � D                          H      
           � D                          P      
           � D                          X      
           � D                          `      
           � D                          h      
           � D                          p      
           � D                          x      
           � D                          �      
           � D                           �                 � D                           �                 � D                             �         
  p      p        p                      � D                             �         
  p      p        p                      � D                           �                 � D                            �                 � D                      !     �                 � D                      "     �                 � D                      #     �                 � D                      $     �                  � D                      %     �   �   !   #INITIALPOWERPARAMS &                 @                   &     '�            #NN '   #AN (   #N_RUN )   #ANT *   #RAT +   #K_0_SCALAR ,   #K_0_TENSOR -   #SCALARPOWERAMP .           �                       '                        �                      (                 
  p      p        p                       �                      )        0         
  p      p        p                       �                      *        X         
  p      p        p                       �                      +        �         
  p      p        p                       �                      ,     �      
            �                      -     �      
            �                      .        �         
  p      p        p                      � D                      /     (   �  "   #REIONIZATIONPARAMS 0              @  @                   0     '(            #REIONIZATION 1   #USE_OPTICAL_DEPTH 2   #REDSHIFT 3   #DELTA_REDSHIFT 4   #FRACTION 5   #OPTICAL_DEPTH 6            � D                      1                        � D                      2                       � D                     3           
            � D                     4           
            � D                     5           
            � D                     6            
           � D                      7           #   #RECOMBINATIONPARAMS 8                 @                   8     '            #RECFAST_FUDGE 9   #RECFAST_FUDGE_HE :   #RECFAST_HESWITCH ;   #RECFAST_HSWITCH <            �                      9            
            �                      :           
            �                       ;                       �                       <                      � D                      =     �    $   #TRANSFERPARAMS >              @  @                   >     '�           #HIGH_PRECISION ?   #NUM_REDSHIFTS @   #KMAX A   #K_PER_LOGINT B   #REDSHIFTS C            � D                      ?                        � D                      @                       � D                     A           
            � D                      B                       � D                     C     �           
  p      p �      p �                    � D  �                   D     
   �     %   
  p      & p    p 
       p 
                     � D                      E        &              � D                      F     0   (  '   #REIONIZATIONHISTORY G              @  @                   G     '0            #TAU_START H   #TAU_COMPLETE I   #AKTHOM J   #FHE K   #WINDOWVARMID L   #WINDOWVARDELTA M            � D                     H            
            � D                     I           
            � D                     J           
            � D                     K           
            � D                     L            
            � D                     M     (      
           � D                      N     X  (              � D                      O     \  )              � D                      P     `  *              � D                     Q     h  +   
           � D                     R     p  ,   
           � D                     S     x  -   
           � D                     T     �  .   
           � D                     U     �  /   
           � D                     V     �  0   
                                    W                                      1                                 X                                      2                                 Y                                      3                               Z                                       [        #     @                          \                   #LENS_CLS%LSAMPLES ]                 @                   ]     '�>           #L0 ^   #L _           �                       ^                        �                       _     �             p      p �      p �            %     @                        `                  
   #BESSI%FLOAT a   #BESSI%ABS b   #BESSI%SQRT c   #BESSI%INT d   #N e   #X f             @                    a     FLOAT           @                    b     ABS           @                    c     SQRT           @                    d     INT       
   @                      e             D @@                      f     
   %     @                      g                  
   #BESSJ0%ABS h   #BESSJ0%SIN i   #BESSJ0%COS j   #BESSJ0%SQRT k   #X l             @                    h     ABS           @                    i     SIN           @                    j     COS           @                    k     SQRT        @@                     l     
      �         fn#fn    �   �   b   uapp(LENSING    L  4   J  PRECISION    �  4   J  MODELPARAMS    �  4   J  AMLUTILS '   �  B     CAMBPARAMS+MODELPARAMS /   *  8   %   CAMBPARAMS%WANTCLS+MODELPARAMS 4   b  8   %   CAMBPARAMS%WANTTRANSFER+MODELPARAMS 3   �  8   %   CAMBPARAMS%WANTSCALARS+MODELPARAMS 3   �  8   %   CAMBPARAMS%WANTTENSORS+MODELPARAMS 3   
  8   %   CAMBPARAMS%WANTVECTORS+MODELPARAMS 1   B  8   %   CAMBPARAMS%DOLENSING+MODELPARAMS 2   z  8   %   CAMBPARAMS%WANT_ZSTAR+MODELPARAMS 2   �  8   %   CAMBPARAMS%WANT_ZDRAG+MODELPARAMS 1   �  8   %   CAMBPARAMS%NONLINEAR+MODELPARAMS -   "  8   %   CAMBPARAMS%MAX_L+MODELPARAMS 4   Z  8   %   CAMBPARAMS%MAX_L_TENSOR+MODELPARAMS 1   �  8   %   CAMBPARAMS%MAX_ETA_K+MODELPARAMS 8   �  8   %   CAMBPARAMS%MAX_ETA_K_TENSOR+MODELPARAMS .     8   %   CAMBPARAMS%OMEGAB+MODELPARAMS .   :  8   %   CAMBPARAMS%OMEGAC+MODELPARAMS .   r  8   %   CAMBPARAMS%OMEGAV+MODELPARAMS .   �  8   %   CAMBPARAMS%OMEGAN+MODELPARAMS *   �  8   %   CAMBPARAMS%H0+MODELPARAMS ,   	  8   %   CAMBPARAMS%TCMB+MODELPARAMS +   R	  8   %   CAMBPARAMS%YHE+MODELPARAMS 7   �	  8   %   CAMBPARAMS%NUM_NU_MASSLESS+MODELPARAMS 6   �	  8   %   CAMBPARAMS%NUM_NU_MASSIVE+MODELPARAMS :   �	  8   %   CAMBPARAMS%NU_MASS_SPLITTINGS+MODELPARAMS ;   2
  8   %   CAMBPARAMS%NU_MASS_EIGENSTATES+MODELPARAMS <   j
  p   %   CAMBPARAMS%NU_MASS_DEGENERACIES+MODELPARAMS 9   �
  p   %   CAMBPARAMS%NU_MASS_FRACTIONS+MODELPARAMS @   J  8   %   CAMBPARAMS%SCALAR_INITIAL_CONDITION+MODELPARAMS ;   �  8   %   CAMBPARAMS%OUTPUTNORMALIZATION+MODELPARAMS <   �  8   %   CAMBPARAMS%ACCURATEPOLARIZATION+MODELPARAMS 2   �  8   %   CAMBPARAMS%ACCURATEBB+MODELPARAMS <   *  8   %   CAMBPARAMS%ACCURATEREIONIZATION+MODELPARAMS 7   b  8   %   CAMBPARAMS%MASSIVENUMETHOD+MODELPARAMS 1   �  P   %   CAMBPARAMS%INITPOWER+MODELPARAMS 0   �  �       INITIALPOWERPARAMS+INITIALPOWER 3   �  8   a   INITIALPOWERPARAMS%NN+INITIALPOWER 3   �  p   a   INITIALPOWERPARAMS%AN+INITIALPOWER 6   /  p   a   INITIALPOWERPARAMS%N_RUN+INITIALPOWER 4   �  p   a   INITIALPOWERPARAMS%ANT+INITIALPOWER 4     p   a   INITIALPOWERPARAMS%RAT+INITIALPOWER ;     8   a   INITIALPOWERPARAMS%K_0_SCALAR+INITIALPOWER ;   �  8   a   INITIALPOWERPARAMS%K_0_TENSOR+INITIALPOWER ?   �  p   a   INITIALPOWERPARAMS%SCALARPOWERAMP+INITIALPOWER -   _  P   %   CAMBPARAMS%REION+MODELPARAMS 0   �  �      REIONIZATIONPARAMS+REIONIZATION =   W  8   %   REIONIZATIONPARAMS%REIONIZATION+REIONIZATION B   �  8   %   REIONIZATIONPARAMS%USE_OPTICAL_DEPTH+REIONIZATION 9   �  8   %   REIONIZATIONPARAMS%REDSHIFT+REIONIZATION ?   �  8   %   REIONIZATIONPARAMS%DELTA_REDSHIFT+REIONIZATION 9   7  8   %   REIONIZATIONPARAMS%FRACTION+REIONIZATION >   o  8   %   REIONIZATIONPARAMS%OPTICAL_DEPTH+REIONIZATION .   �  Q   %   CAMBPARAMS%RECOMB+MODELPARAMS 2   �  �       RECOMBINATIONPARAMS+RECOMBINATION @   �  8   a   RECOMBINATIONPARAMS%RECFAST_FUDGE+RECOMBINATION C   �  8   a   RECOMBINATIONPARAMS%RECFAST_FUDGE_HE+RECOMBINATION C   �  8   a   RECOMBINATIONPARAMS%RECFAST_HESWITCH+RECOMBINATION B   0  8   a   RECOMBINATIONPARAMS%RECFAST_HSWITCH+RECOMBINATION 0   h  L   %   CAMBPARAMS%TRANSFER+MODELPARAMS +   �  �      TRANSFERPARAMS+MODELPARAMS :   B  8   %   TRANSFERPARAMS%HIGH_PRECISION+MODELPARAMS 9   z  8   %   TRANSFERPARAMS%NUM_REDSHIFTS+MODELPARAMS 0   �  8   %   TRANSFERPARAMS%KMAX+MODELPARAMS 8   �  8   %   TRANSFERPARAMS%K_PER_LOGINT+MODELPARAMS 5   "  p   %   TRANSFERPARAMS%REDSHIFTS+MODELPARAMS >   �  |   %   CAMBPARAMS%INITIALCONDITIONVECTOR+MODELPARAMS 5     8   %   CAMBPARAMS%ONLYTRANSFERS+MODELPARAMS 1   F  Q   %   CAMBPARAMS%REIONHIST+MODELPARAMS 1   �  �      REIONIZATIONHISTORY+REIONIZATION ;   /  8   %   REIONIZATIONHISTORY%TAU_START+REIONIZATION >   g  8   %   REIONIZATIONHISTORY%TAU_COMPLETE+REIONIZATION 8   �  8   %   REIONIZATIONHISTORY%AKTHOM+REIONIZATION 5   �  8   %   REIONIZATIONHISTORY%FHE+REIONIZATION >     8   %   REIONIZATIONHISTORY%WINDOWVARMID+REIONIZATION @   G  8   %   REIONIZATIONHISTORY%WINDOWVARDELTA+REIONIZATION ,     8   %   CAMBPARAMS%FLAT+MODELPARAMS .   �  8   %   CAMBPARAMS%CLOSED+MODELPARAMS ,   �  8   %   CAMBPARAMS%OPEN+MODELPARAMS .   '  8   %   CAMBPARAMS%OMEGAK+MODELPARAMS ,   _  8   %   CAMBPARAMS%CURV+MODELPARAMS )   �  8   %   CAMBPARAMS%R+MODELPARAMS -   �  8   %   CAMBPARAMS%KSIGN+MODELPARAMS ,     8   %   CAMBPARAMS%TAU0+MODELPARAMS ,   ?  8   %   CAMBPARAMS%CHI0+MODELPARAMS )   w  U       LENSING_METHOD_CURV_CORR )   �  U       LENSING_METHOD_FLAT_CORR (   !  U       LENSING_METHOD_HARMONIC    v  0       LENSING_METHOD )   �  0       LENSING_INCLUDES_TENSORS    �  S       LENS_CLS *   )  K      LENS_CLS%LSAMPLES+LVALUES -   t  8   a   LENS_CLS%LSAMPLES%L0+LVALUES ,   �  p   a   LENS_CLS%LSAMPLES%L+LVALUES      �       BESSI    �  2      BESSI%FLOAT    �  0      BESSI%ABS      1      BESSI%SQRT    <  0      BESSI%INT    l  0   a   BESSI%N    �  0   a   BESSI%X    �  �       BESSJ0    T   0      BESSJ0%ABS    �   0      BESSJ0%SIN    �   0      BESSJ0%COS    �   1      BESSJ0%SQRT    !  0   a   BESSJ0%X 