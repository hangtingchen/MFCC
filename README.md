# MFCC

C/C++ code to extract MFCC or FBank features from wav files.

# masterCPLus should be used. The mater branch may not be updated in time.

# Install

1. Download following code from my GitHub and put these three directories under the same root directory.
    1. https://github.com/hangtingchen/BasicAudioToolBox
    2. https://github.com/hangtingchen/MFCC
    3. https://github.com/hangtingchen/inih
2. Switch the branch from master to masterCPlus in BasicAudioToolBox and MFCC. (The master branch is okay but it is highly recommended to use masterCPlus)
3. For Linux, use `make` in `MFCC` directory to construct the program. The target is `mfcc`.
4. For Windows, 
    1. Add these files to your visual studio.
    2. Add path of `BasicAudioToolBox` and `inih` to VC++ catalog.
    3. Add `_CRT_SECURE_NO_WARNINGS` to the preprocessor definition.
    4. Enable `openmp`.(if you don't want to use multiple threads, ignore this step)
    5. Generate the exe.

# Example

Run an example first. Enter `MFCC` directory.
```shell
mfcc example/config.ini
```

The screen will display followings,

```
finish reading the config.ini
FFT passband 2 to 1024 out of 1 to 1024
Mel passband 15.986084 to 3923.357581
Thread ID 1
Convert example/a001_0_30.wav to example/a001_0_30.fbank
including : 
MFCCNum 40
energyFlag 1
zeroCrossingFlag 1
brightFlag 1
subBandEFlag 8
the frame feature dimension is 408
Sample Rate : 44100
Number of channels : 2
Each sample's size in byte : 3
Each container's size in byte : 3
Number of samples : 1323001
total coef size: 612000
post-processing...
writing the doc...
```

These message indicts following steps to extract MFCC,

1. The program will first read `config.ini`. 
2. The Mel filter bands will be generated according to the setting in config.
3. Read WAV file.
4. Count and extract the features.
5. Write the feature file.

# Usage

You should prepare WAV files first(.wav), then a config file and a list. The following introduces the specific format.

## WAV

The WAV should be PCM encoding and has a standard 44-byte head. (It's okay if there exists additional information chunk between 44-byte head and data chunk.) Don't worry about this, for most of WAV files are satisfied. If the wav is transformed from mp3 or other format, have a look at the head first. 

You should check WAV if No new information displayed in the screen for a long time while running the program. 

## config

The config is read using `inih`(https://github.com/benhoyt/inih). A standard config is `example/config.ini` under `MFCC` directory. The options listed should all be included. No promise for the program if you skip some options.

Website(https://github.com/hangtingchen/inih/blob/master/examples/test.ini) gives details about `inih`.

### `Frame` section

This section sets options related to pre-process.

| Key | Value | Notes |
| :------| :------ | :------ |
| sampleRate | 8000/16000/44100/others | The sample rate should be set at first, which means the wav files should have same sample rate. |
| lowpassfre | >0 && <hipassfre && <sampleRate | The min frequency |
| hipassfre | >0 && >lowpassfre && <=sampleRate | The max frequency |
| preemphasise | =0 no preemphasise; 0-1 preemphasise the signal | The coefficient of preemphasise |
| zeroMeanSigFlag | 0/1 | whether to make input signal have a zero mean |
| wlen | wlen=(wlenInTime(ms))*sampleRate/1000 | the number of samples of window |
| inc | inc=(incInTime(ms))*sampleRate/1000 | the number of samples of window shift |
| vecNum | =1(mono/double channel WAV);=2(double channel WAV);=4(double channel WAV) | the channels of the output feature |

### `MFCC` section 

This section controls the Mel filter settings.

| Key | Value | Notes |
| :------| :------ | :------ |
| fbankFlag | 0/1 | Extract MFCC or fbank |
| bankNum | >0 | The number of Mel filters |
| MFCCNum | >0 && <=bankNum | The number of MFCC; If fbankFlag=1, this option has no effect |
| MFCC0thFlag | 0/1 | Whether to include MFCC0th; If fbankFlag=1, this option has no effect |

### `Others` section

This section controls some other features.

| Key | Value | Notes |
| :------| :------ | :------ |
| energyFlag | 0/1 | Whether to include average energy |
| zeroCrossingFlag | 0/1 | Whether to include average zero crossing rate |
| brightFlag | 0/1 | Whether to include brightness |
| subBandEFlag | =0, No subBand energy; >0 set the number of subband | |
| fftLength |  >=0 | Output fft. This is only for debugging. |

### `Regression` section

This section controls post-process.

| Key | Value | Notes |
| :------| :------ | :------ |
| znormFlag | 0/1 | Whether to do z-norm in each dimension within the single audio |
| regreOrder | =1,no diff;=2 first order diff; and so on | (The degree of diff) + 1 |
| delwin | >0 | The context length of diff; If regreOrder=1, this option has no effect |

### `IO` section

This section controls how to read file and store features.

| Key | Value | Notes |
| :------| :------ | :------ |
| fileList | | The position of list of files |
| saveType | f/e/n/b | csv(double)/csv(scientific)/npy(numpy)/binary |
| numThreads | >0 | The number of threads.(Only when openmp is supported.) |

## file list

Each line in the list should include the source WAV and target feature, separated by `tab`. Please refer to `example/fileList.txt`.

# Notes

+ The main interface is `MFCCWapper` in `mfcc.c` in case you want to extract feature in your own code or you don't want to control so many options. 
+ The program has set up `hmath` and `hsigProcess` according to HTK.
+ Recommend to use `masterCPlus` branch.
+ The anomaly detection is not perfect. Please pay attention to your config and file list.
+ Though the code in `masterCPlus` is C++, the style is still C,and so is the memory application and release.