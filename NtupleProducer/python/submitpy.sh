
#!/bin/bash

#python submit.py --sample dy_140 --tkptcut 2.0 --submit
#python submit.py --sample tt_140 --tkptcut 2.0 --submit
#python submit.py --sample qcdmu_140 --tkptcut 2.0 --submit
#python submit.py --sample mugun_140 --tkptcut 2.0 --metRate --submit

python submit.py --sample dy_140 --tkptcut 2.5 --submit
python submit.py --sample tt_140 --tkptcut 2.5 --submit
python submit.py --sample qcdmu_140 --tkptcut 2.5 --submit
python submit.py --sample mugun_140 --tkptcut 2.5 --metRate --submit

python submit.py --sample dy_140 --tkptcut 3.0 --submit
python submit.py --sample tt_140 --tkptcut 3.0 --submit
python submit.py --sample qcdmu_140 --tkptcut 3.0 --submit
python submit.py --sample mugun_140 --tkptcut 3.0 --metRate --submit

python submit.py --sample dy_140 --tkptcut 0.0 --submit
python submit.py --sample tt_140 --tkptcut 0.0 --submit
python submit.py --sample qcdmu_140 --tkptcut 0.0 --submit
python submit.py --sample mugun_140 --tkptcut 0.0 --metRate --submit