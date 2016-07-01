#!/bin/bash

if  [ ! -f profile.out ]; then
  python -m cProfile -o profile.out main_1.py ../out.xml
fi
python -c "import pstats; p=pstats.Stats('profile.out'); p.sort_stats('time').print_stats()" 
