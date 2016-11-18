#!/bin/bash

ALL="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 \
26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48   \
49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71   \
72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94   \
95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 \
114 115 116 117 118 119 201 202 203 204 205 206 207 208 209 210 211    \
212 213 214 215 216 217 218 219 220 221 222 223 224 225 226 227 228    \
229 230 231 232 233 234 235 236 237 238 239 240 241 242 243 244 245    \
246 247 248 249 250 251 252 253 254 255 256 257 258 259 260 261 262    \
263 264 265 266 267 268 269 270 271 272 273 274 275 276 277 278 279    \
280 281 282 283 284 285 286 287 288 289 290 291 292 293 294 295 296    \
297 298 299 300 301 302 303 304 305 306 307 308 309 310 311 312 313    \
314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330    \
331 332 333 334 335 336 337 338 339 340 341 342 343 344 345 346 347    \
348 349 350 351 352 353 354 355 356 357 358 359 360 361 362 363 364    \
365 366 367 368 369 370 371 372 373 374 375 376 377 378 379 380 381    \
382 383 384 385 386 387 388 389 390 391 392 393 394 395"

EQ="6 7 8 9 14 26 27 28 32 39 40 41 42 46 47 48 49 50 51 52 53 54  \
55 56 60 61 62 63 68 69 71 73 74 75 77 78 79 80 81 87 99 107 109 111 \
112 114 119 216 217 219 235 248 252 254 262 263 265 269 316 317 318  \
319 320 321 322 325 335 336 338 344 345 348 353 355 367 373 375 376  \
377 378 381 382 383 390 393 394 395"

INEQ="10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 29 30 31 32 33    \
34 35 36 37 43 44 57 58 59 64 65 66 67 70 71 72 73 74 75 76 83 84 85 86 88 \
89 90 91 92 93 95 96 97 98 100 101 102 103 104 105 106 108 109 113   \
114 116 117 118 215 217 218 220 221 222 223 224 225 226 227 228 230  \
231 232 233 234 236 237 238 239 248 249 250 251 253 262 263 264 268  \
270 277 278 279 280 284 285 315 323 324 325 326 327 329 330 331 332  \
337 339 340 341 342 343 346 347 349 353 354 356 359 360 361 362 363  \
364 365 366 367 369 372 374 376 380 381 382 384 385 386 387 388 389  \
390 392 393"

CONSTRAINED="6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24    \
26 27 28 29 30 31 32 33 34 35 36 37 39 40 41 42 43 44 46 47 48 49 50 \
51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 \
75 76 77 78 79 80 81 83 84 85 86 87 88 89 90 91 92 93 95 96 97 98 99    \
100 101 102 103 104 105 106 107 108 109 111 112 113 114 116 117 118  \
119 215 216 217 218 219 220 221 222 223 224 225 226 227 228 230 231  \
232 233 234 235 236 237 238 239 248 249 250 251 252 253 254 262 263  \
264 265 268 269 270 277 278 279 280 284 285 315 316 317 318 319 320  \
321 322 323 324 325 326 327 329 330 331 332 335 336 337 338 339 340  \
341 342 343 344 345 346 347 348 349 353 354 355 356 359 360 361 362  \
363 364 365 366 367 369 372 373 374 375 376 377 378 380 381 382 383  \
384 385 386 387 388 389 390 392 393 394 395"

CONSTRDER="6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24    \
26 27 28 29 30 31 32 33 34 35 36 37 39 40 41 42 43 44 46 47 48 49 50 \
51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 \
75 76 77 78 79 80 81 83 84 85 86 87 88 89 90 91 92 93 95 96 97 98 99    \
100 101 102 103 104 105 106 107 108 109 111 112 113 114 116 117 118  \
119 215 216 217 218 219 220 221 222 223 224 225 226 227 228 230 231  \
232 233 234 235 236 237 238 239 248 249 250 251 252 253 254 262 263  \
264 265 268 269 270 277 278 279 280 284 285 315 316 317 318 319 320  \
321 322 323 324 325 326 327 329 330 331 335 336 337 338 339 340  \
341 342 343 344 345 346 347 348 349 353 354 355 359 360 361   \
367 372 373 374 375 376 377 378 380 381 382 383  \
384 385 386 387 388 389 394 395"

NOCONSTRDER="332 356 362 363 364 365 366 369 390 392 393"

if (test $# = 1 ) then
    
    EXECNAME=$1;
    
    ulimit -St 1800 # 30 min of CPU time per problem

    rm -f analysis-all
    
    for i in $CONSTRAINED; do
        echo $i | ./$EXECNAME

        cat runhs.out >> analysis-all
        rm -f runhs.out

        # Only for ALGENCAN's tests
        #
        #awk '{ printf "%15d %9.2f\n",$20,$25 }' algencan-tabline.out >> \
        #    analysis-all

    done;

else

    echo -e "\n\nUSO: ./runhstests.sh <NOMEEXEC>\n\n";

fi
