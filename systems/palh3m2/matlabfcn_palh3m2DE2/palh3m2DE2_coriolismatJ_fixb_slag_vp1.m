% Calculate matrix of centrifugal and coriolis load on the joints for
% palh3m2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [18x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,BC,BE,BG,DC,DT2,EP,GH,GP,HW,OT1,T1A,T1T2,phi1,phi2,phi410,phi78,phi79]';
% m [9x1]
%   mass of all robot links (including the base)
% rSges [9x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [9x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 04:24
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = palh3m2DE2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(18,1),zeros(9,1),zeros(9,3),zeros(9,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh3m2DE2_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'palh3m2DE2_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [18 1]), ...
  'palh3m2DE2_coriolismatJ_fixb_slag_vp1: pkin has to be [18x1] (double)');
assert(isreal(m) && all(size(m) == [9 1]), ...
  'palh3m2DE2_coriolismatJ_fixb_slag_vp1: m has to be [9x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [9,3]), ...
  'palh3m2DE2_coriolismatJ_fixb_slag_vp1: rSges has to be [9x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [9 6]), ...
  'palh3m2DE2_coriolismatJ_fixb_slag_vp1: Icges has to be [9x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 04:15:24
% EndTime: 2020-05-07 04:15:48
% DurationCPUTime: 24.92s
% Computational Cost: add. (72241->595), mult. (26962->778), div. (4974->2), fcn. (29206->133), ass. (0->465)
t3188 = sin(qJ(2));
t3192 = cos(qJ(2));
t3428 = t3192 * t3188;
t3381 = rSges(9,2) * t3428;
t3322 = rSges(9,1) * t3381;
t3076 = -0.2e1 * t3322;
t3478 = -qJ(2) - pkin(15);
t3155 = pkin(18) - t3478;
t3110 = pkin(17) + qJ(3) + t3155;
t3306 = atan2(-sin(t3110), -cos(t3110));
t3117 = sin(t3155);
t3118 = cos(t3155);
t3021 = atan2(t3117, -t3118);
t3426 = pkin(17) - t3021;
t2912 = -t3306 + t3426;
t2907 = sin(t2912);
t2908 = cos(t2912);
t3455 = t2907 * t2908;
t3201 = rSges(9,2) ^ 2;
t3203 = rSges(9,1) ^ 2;
t3167 = -t3201 + t3203;
t3174 = t3192 ^ 2;
t3107 = t3167 * t3174;
t3137 = t3203 / 0.2e1 - t3201 / 0.2e1;
t3613 = t3137 - t3107;
t3367 = (t3076 - t3613) * t3455;
t3315 = m(9) * t3367;
t3640 = 0.2e1 * t3315;
t3184 = sin(pkin(16));
t3185 = cos(pkin(16));
t3189 = sin(pkin(15));
t3194 = cos(pkin(15));
t3027 = t3184 * t3194 + t3185 * t3189;
t3028 = -t3184 * t3189 + t3185 * t3194;
t3187 = sin(qJ(3));
t3191 = cos(qJ(3));
t2915 = -t3027 * t3191 - t3187 * t3028;
t3267 = t3187 * t3027 - t3028 * t3191;
t2875 = t2915 * t3192 + t3188 * t3267;
t3170 = pkin(17) + pkin(18);
t3139 = sin(t3170);
t3140 = cos(t3170);
t3477 = qJ(2) + qJ(3);
t3611 = -t3188 * t2915 + t3267 * t3192;
t2767 = sin(atan2(t2875 * t3140 + t3139 * t3611, t2875 * t3139 - t3611 * t3140) + t3477);
t3263 = pkin(16) + t3110;
t2980 = atan2(-sin(t3263), cos(t3263));
t3474 = qJ(3) + t2980;
t2977 = qJ(2) + t3474;
t2971 = qJ(1) + t2977;
t2972 = qJ(1) - t2977;
t3193 = cos(qJ(1));
t3479 = t3193 * rSges(6,1);
t3541 = sin(qJ(1));
t3161 = t3541 * rSges(6,2);
t3526 = t3161 / 0.2e1;
t3044 = t3526 + t3479 / 0.2e1;
t3045 = t3526 - t3479 / 0.2e1;
t3162 = t3193 * rSges(6,2);
t3394 = t3541 * rSges(6,1);
t3314 = -t3394 / 0.2e1;
t3046 = t3314 + t3162 / 0.2e1;
t3047 = t3314 - t3162 / 0.2e1;
t3061 = -t3161 + t3479;
t3062 = -t3161 - t3479;
t3063 = t3162 + t3394;
t3064 = t3162 - t3394;
t3177 = qJ(1) + qJ(4);
t3144 = sin(t3177);
t3178 = qJ(1) - qJ(4);
t3145 = sin(t3178);
t3151 = cos(t3177);
t3152 = cos(t3178);
t3473 = t2980 - qJ(4);
t2976 = qJ(3) + t3473;
t2970 = qJ(2) + t2976;
t2964 = qJ(1) + t2970;
t2950 = cos(t2964);
t3472 = t2980 + qJ(4);
t2975 = qJ(3) + t3472;
t2969 = qJ(2) + t2975;
t2965 = qJ(1) - t2969;
t2951 = cos(t2965);
t3418 = t2950 + t2951;
t2946 = sin(t2964);
t2947 = sin(t2965);
t3419 = t2946 + t2947;
t2962 = qJ(1) + t2969;
t2944 = sin(t2962);
t2963 = qJ(1) - t2970;
t2945 = sin(t2963);
t3420 = t2944 + t2945;
t2948 = cos(t2962);
t2949 = cos(t2963);
t3621 = t2948 + t2949;
t2764 = -t3061 * t3145 - t3144 * t3062 + t3063 * t3152 + t3064 * t3151 + t3418 * t3047 + t3621 * t3046 - t3419 * t3045 + t3420 * t3044 + ((-cos(t2971) + cos(t2972)) * t3193 + (-sin(t2971) + sin(t2972)) * t3541) * rSges(6,3);
t3143 = sin(t3477);
t3506 = pkin(4) * t3143;
t3390 = t2764 * t3506;
t3208 = t2767 * t3390;
t3634 = -m(6) / 0.4e1;
t3639 = rSges(6,2) * t3634;
t3439 = t3117 ^ 2 / t3118 ^ 2;
t3338 = 0.1e1 + t3439;
t3019 = 0.1e1 / t3338;
t3616 = t3019 * t3338;
t3631 = t3616 - 0.1e1;
t2861 = t3631 - 0.1e1;
t2909 = -qJ(2) + t2912;
t2903 = sin(t2909);
t2904 = cos(t2909);
t3545 = m(9) * rSges(9,3);
t3628 = (rSges(9,2) * t3545 - Icges(9,6)) * t2903 / 0.2e1 - (-rSges(9,1) * t3545 + Icges(9,5)) * t2904 / 0.2e1;
t3637 = t2861 * t3628;
t3635 = pkin(4) / 0.2e1;
t3633 = m(6) / 0.4e1;
t3632 = cos(t2969) / 0.2e1;
t3559 = pkin(1) * m(6);
t3221 = t3559 / 0.4e1;
t3557 = pkin(4) * m(6);
t3405 = t3557 / 0.4e1;
t3630 = rSges(6,1) * t3405;
t3175 = qJ(4) + qJ(2);
t3148 = cos(t3175);
t3176 = -qJ(4) + qJ(2);
t3149 = cos(t3176);
t3156 = qJ(3) + t3175;
t3127 = cos(t3156);
t3157 = qJ(3) + t3176;
t3128 = cos(t3157);
t3542 = rSges(6,2) * m(6);
t3397 = t3542 / 0.4e1;
t3347 = pkin(4) * t3397;
t3399 = -t3542 / 0.4e1;
t3349 = pkin(4) * t3399;
t3122 = sin(t3156);
t3123 = sin(t3157);
t3614 = t3122 + t3123;
t3291 = t3127 * t3347 + t3128 * t3349 + t3614 * t3630;
t3400 = rSges(6,1) * t3634;
t3141 = sin(t3175);
t3142 = sin(t3176);
t3620 = t3141 + t3142;
t3629 = t3291 + (t3148 * t3399 + t3149 * t3397 + t3400 * t3620) * pkin(1);
t3626 = t2945 / 0.4e1 + t2944 / 0.4e1 + t3144 / 0.2e1;
t3624 = -0.2e1 * t3628;
t3555 = -rSges(9,1) / 0.2e1;
t3432 = t3187 * t3191;
t3205 = pkin(4) ^ 2;
t3492 = t3205 * m(5);
t3101 = t3432 * t3492;
t3540 = rSges(9,1) * rSges(9,2);
t3407 = m(9) * t3540;
t3037 = t3101 - t3407;
t3172 = t3188 ^ 2;
t3595 = t3172 - t3174;
t3622 = t3595 * t3037;
t3034 = t3188 * t3187 - t3192 * t3191;
t3035 = t3192 * t3187 + t3188 * t3191;
t2929 = t3034 * t3194 + t3189 * t3035;
t3594 = -t3189 * t3034 + t3035 * t3194;
t2878 = t2929 * t3185 + t3594 * t3184;
t2902 = t2904 ^ 2;
t3401 = m(9) * rSges(9,2);
t2854 = rSges(9,1) * t2902 * t3401;
t3433 = t3167 * t2904;
t3362 = t2903 * t3433;
t3316 = m(9) * t3362;
t3615 = -t2854 - t3101 + t3316;
t2906 = t2908 ^ 2;
t3361 = t3167 * t3428;
t3164 = rSges(9,2) * t3174;
t3138 = rSges(9,1) * t3164;
t3598 = 0.2e1 * t3138 - t3540;
t3255 = t3361 + t3598;
t3165 = rSges(9,1) * t3174;
t3513 = pkin(1) * (rSges(9,1) - t3165 + t3381);
t3457 = (-t2907 * t3255 - t3513) * t2907;
t3612 = -t3255 * t2906 - t3457;
t3283 = t2949 / 0.4e1 + t2948 / 0.4e1 + t3151 / 0.2e1;
t3159 = qJ(1) + t3477;
t3125 = sin(t3159);
t3509 = pkin(4) * t3125;
t3077 = t3509 / 0.2e1;
t3160 = qJ(1) - t3477;
t3126 = sin(t3160);
t3508 = pkin(4) * t3126;
t3572 = -t3508 + 0.2e1 * t3077;
t3130 = cos(t3160);
t3507 = pkin(4) * t3130;
t3609 = t3507 / 0.4e1;
t3568 = -pkin(1) / 0.2e1;
t3606 = m(8) * t3568;
t3271 = (t2907 ^ 2 - t2906) * (0.2e1 * t3137 * t3428 + t3598);
t2901 = t2903 ^ 2;
t3601 = (t2901 - t2902) * t3540;
t3150 = cos(t3477);
t3488 = rSges(4,2) * t3192;
t3441 = (-rSges(4,1) * t3188 + t3488) * t3150;
t3560 = pkin(1) * m(4);
t3007 = t3441 * t3560;
t3476 = pkin(1) * pkin(4) * m(5);
t3120 = t3187 * t3476;
t3112 = -0.4e1 * t3120;
t3119 = t3191 * t3476;
t3173 = t3191 ^ 2;
t3411 = t3187 ^ 2 - t3173;
t3289 = t3411 * t3492;
t3383 = rSges(9,1) * t3428;
t3498 = m(9) * t2908;
t3333 = pkin(1) * (-rSges(9,2) + t3164 + t3383) * t3498;
t3412 = t3120 - t3333;
t3525 = t3174 / 0.4e1;
t3600 = -(t3119 + t3289) * t3428 + (0.8e1 * t3101 + t3112) * t3525 + t3007 + t3412 + t3615;
t3129 = cos(t3159);
t3179 = qJ(1) + qJ(2);
t3146 = sin(t3179);
t3153 = cos(t3179);
t3378 = pkin(1) * t3405;
t3379 = -pkin(4) * t3559 / 0.4e1;
t3599 = t3153 * t3125 * t3378 + t3146 * t3129 * t3379;
t3195 = pkin(10) + rSges(6,3);
t3163 = t3195 * rSges(6,1);
t3543 = rSges(6,2) * pkin(8);
t3597 = t3163 - t3543;
t3596 = t3163 + t3543;
t3495 = t3167 * m(9);
t3102 = 0.4e1 * t3495;
t3578 = 0.2e1 * t2909;
t3591 = sin(t3578) * (t3102 - (4 * Icges(9,1)) + (4 * Icges(9,2))) / 0.8e1 - (-Icges(9,4) + t3407) * cos(t3578);
t3486 = rSges(9,2) * t2903;
t3499 = m(9) * t2907;
t3510 = pkin(1) * t3188;
t3590 = t3631 * t3333 + pkin(1) * (-rSges(9,1) * t3172 + t3165 - 0.2e1 * t3381) * t3499 + m(9) * t3486 * t3510;
t3282 = t2951 / 0.4e1 + t2950 / 0.4e1 - t3152 / 0.2e1;
t3284 = t2947 / 0.4e1 - t3145 / 0.2e1 + t2946 / 0.4e1;
t3197 = m(5) + m(6);
t3589 = pkin(4) * t3197 + m(4) * rSges(4,1);
t3586 = -t3284 * rSges(6,1) + t3282 * rSges(6,2);
t3585 = pkin(1) ^ 2;
t3584 = -4 * Icges(6,5);
t3583 = -0.2e1 * qJ(2);
t3582 = 0.2e1 * qJ(3);
t3580 = -0.2e1 * (-t3626 - t3284) * rSges(6,2) - 0.2e1 * (-t3282 + t3283) * rSges(6,1);
t3259 = t3626 * rSges(6,1) + t3283 * rSges(6,2);
t3579 = 0.2e1 * t3259 + 0.2e1 * t3586;
t3575 = 0.2e1 * t2969;
t3574 = 0.2e1 * t2970;
t3008 = qJ(2) + t3021;
t3570 = 0.2e1 * t3008;
t3158 = pkin(14) + t3478;
t3569 = 0.2e1 * t3158;
t3567 = pkin(1) / 0.2e1;
t3564 = -m(6) / 0.2e1;
t3562 = m(6) / 0.2e1;
t3558 = pkin(3) * m(9);
t3556 = pkin(8) * m(6);
t3554 = pkin(4) * rSges(6,1);
t3553 = pkin(4) * rSges(6,2);
t3551 = m(4) * rSges(4,2);
t3550 = m(4) * rSges(4,3);
t3548 = m(5) * rSges(5,3);
t3547 = m(7) * rSges(7,2);
t3546 = m(7) * rSges(7,3);
t3544 = rSges(3,2) * m(3);
t3109 = pkin(4) * t3129;
t3530 = -t3109 / 0.2e1;
t3199 = 0.2e1 * qJ(2);
t3417 = t2980 + t3199;
t3375 = qJ(3) + t3417;
t2967 = qJ(4) + t3375;
t3524 = cos(t2967);
t2968 = -qJ(4) + t3375;
t3523 = cos(t2968);
t3522 = cos(t2975);
t3521 = cos(t2976);
t3520 = cos(t3158);
t3519 = sin(t2967);
t3518 = sin(t2968);
t3517 = sin(t2975);
t3516 = sin(t2976);
t3484 = rSges(9,2) * t3172;
t3514 = pkin(1) * (-t3164 - 0.2e1 * t3383 + t3484);
t3180 = qJ(1) - qJ(2);
t3147 = sin(t3180);
t3512 = pkin(1) * t3147;
t3511 = pkin(1) * t3153;
t3505 = pkin(4) * t3150;
t3504 = m(6) * qJD(2);
t3501 = m(6) * t3205;
t3088 = t3123 * t3554;
t3303 = t3128 * t3553 - t3088;
t3087 = t3122 * t3554;
t3413 = t3127 * t3553 + t3087;
t2914 = -t3303 + t3413;
t3497 = (((-t3148 + t3149) * rSges(6,2) - t3620 * rSges(6,1)) * pkin(1) + t2914) * m(6);
t3496 = t2914 * m(6);
t3494 = (rSges(6,1) ^ 2 - rSges(6,2) ^ 2) * m(6);
t3169 = t3192 * pkin(1);
t3491 = rSges(9,1) * t2903;
t3490 = rSges(9,1) * t2904;
t3489 = rSges(4,2) * t3188;
t3487 = rSges(6,2) * t3195;
t3485 = rSges(9,2) * t2904;
t3482 = pkin(12) * t3631;
t3481 = sin(qJ(4)) * rSges(6,1);
t3480 = cos(qJ(4)) * rSges(6,2);
t3471 = t3199 + qJ(3);
t3470 = 0.4e1 * m(6);
t3321 = t2929 * t3184 - t3185 * t3594;
t2770 = sin(atan2(t2878 * t3139 + t3140 * t3321, -t2878 * t3140 + t3139 * t3321) + t3477);
t3244 = t2770 * t3390;
t2736 = t3208 * t3634 + t3244 * t3633;
t3468 = t2736 * qJD(4);
t3262 = pkin(1) * t3589;
t3023 = sin(t3471) * t3262;
t3427 = pkin(1) * t3551;
t3092 = cos(t3471) * t3427;
t3124 = sin(t3158);
t3216 = rSges(6,2) * t3633;
t3374 = t3582 + t3417;
t3307 = qJ(4) + t3374;
t3272 = sin(t3307);
t3308 = -qJ(4) + t3374;
t3273 = sin(t3308);
t3274 = cos(t3307);
t3275 = cos(t3308);
t3290 = -m(5) * rSges(5,2) + t3195 * m(6);
t3385 = -m(9) / 0.2e1;
t3323 = rSges(9,2) * t3385;
t3293 = pkin(3) * t3323;
t3311 = m(5) * rSges(5,1) + t3556;
t3319 = 0.2e1 * t3477;
t3368 = sin(t3472);
t3369 = -sin(t3473);
t3371 = cos(t3472);
t3372 = cos(t3473);
t3207 = ((rSges(4,1) ^ 2 - rSges(4,2) ^ 2) * m(4) - Icges(4,1) + Icges(4,2) + t3197 * t3205) * sin(t3319) / 0.2e1 + cos(t3306) * t3293 - (-rSges(4,1) * t3551 + Icges(4,4)) * cos(t3319) + pkin(3) * rSges(9,1) * sin(t3306) * t3385 + (-sin(t2980) + sin(t3374)) * t3311 * t3635 + t3372 * t3405 * rSges(6,2) + (t3368 * t3634 + (t3272 + t3273) * t3633) * t3554 + t3369 * t3630 + (-t3589 * t3143 - t3150 * t3551) * pkin(12) + (-(-cos(t2980) + cos(t3374)) * t3290 / 0.2e1 + t3274 * t3216 + (t3275 + t3371) * t3639) * pkin(4);
t3231 = (-0.3e1 + (0.2e1 + 0.2e1 * t3439) * t3019) * t3558;
t3233 = pkin(1) * m(9) * (t3631 - 0.2e1);
t3247 = t3583 + t2912;
t3240 = sin(t3247);
t3248 = t3583 - t3306 - 0.2e1 * t3021 + 0.2e1 * pkin(17);
t3241 = sin(t3248);
t3242 = cos(t3247);
t3243 = cos(t3248);
t3266 = t3616 - 0.2e1;
t3245 = t3266 * t3606;
t3252 = t3616 * t3606;
t3402 = m(9) * t3491;
t3328 = t2861 * t3402;
t3337 = m(4) + m(8) + m(9) + t3197;
t3377 = -qJ(2) + t3426;
t3382 = t2861 * t3485;
t3462 = t3631 * t2908;
t3391 = m(9) * t3462;
t3392 = t3631 * t3499;
t3403 = m(8) * t3482;
t3409 = t3199 + t3021;
t2739 = -t3092 + ((rSges(3,1) ^ 2 - rSges(3,2) ^ 2) * m(3) - Icges(3,1) + Icges(3,2) + t3337 * t3585) * sin(t3199) / 0.2e1 + (rSges(7,1) * t3547 - Icges(7,4)) * cos(t3569) - (-rSges(3,1) * t3544 + Icges(3,4)) * cos(t3199) - ((rSges(7,1) ^ 2 - rSges(7,2) ^ 2) * m(7) - Icges(7,1) + Icges(7,2)) * sin(t3569) / 0.2e1 + t3207 - t3023 + rSges(6,2) * t3522 * t3221 + (sin(t3375) - sin(t3474)) * t3311 * t3568 + (cos(t3375) - cos(t3474)) * t3290 * t3567 + (t3241 * t3231 + t3240 * t3233) * t3555 + t3591 * t2861 + (m(7) * rSges(7,1) * t3124 - t3520 * t3547) * pkin(6) + (-cos(t3008) * t3403 + cos(t3021) * t3252 + cos(t3409) * t3245) * rSges(8,2) + (-sin(t3008) * t3403 + sin(t3021) * t3252 + sin(t3409) * t3245) * rSges(8,1) + ((-m(8) * rSges(8,1) * rSges(8,2) + Icges(8,4)) * cos(t3570) - ((rSges(8,1) ^ 2 - rSges(8,2) ^ 2) * m(8) - Icges(8,1) + Icges(8,2)) * sin(t3570) / 0.2e1 + pkin(3) ^ 2 * m(9) * sin(0.2e1 * t3377) / 0.2e1) * t3631 + (t3391 * t3567 + t3242 * t3233 / 0.2e1 + t3243 * t3231 / 0.2e1) * rSges(9,2) + (t3523 * t3216 + t3392 * t3555 + (t3521 + t3524) * t3639) * pkin(1) + (-t3328 + (t3337 * pkin(1) + rSges(3,1) * m(3)) * t3188 + t3192 * t3544 + m(9) * t3382) * pkin(12) + (sin(t3377) * t3482 + (sin(t3426) * t3616 + t3266 * sin(t3583 + t3426)) * t3567) * t3558 + ((t3518 + t3519) * pkin(1) * t3634 + (t3516 + t3517) * t3221) * rSges(6,1);
t3467 = t2739 * qJD(1);
t3324 = m(9) * t3555;
t3326 = t2904 * t3401;
t2740 = -t3092 / 0.2e1 + t3207 - t3191 * t3427 / 0.2e1 - pkin(3) * t3241 * t3324 - t3023 / 0.2e1 + t3243 * t3293 - t3187 * t3262 / 0.2e1 - t3591 + (-t3326 + t3402) * pkin(12) + (-(t2907 + t3240) * t3324 + (t3242 + t2908) * t3323) * pkin(1);
t3466 = t2740 * qJD(1);
t2954 = sin(t2969);
t2955 = sin(t2970);
t2959 = cos(t2970);
t3196 = rSges(6,1) * pkin(8);
t3040 = (t3196 + t3487) * m(6) - Icges(6,6);
t3041 = (t3196 - t3487) * m(6) + Icges(6,6);
t3198 = 0.2e1 * qJ(4);
t3376 = 0.2e1 * t2980 + t3199 + t3582;
t3309 = qJ(4) + t3376;
t3310 = -qJ(4) + t3376;
t3320 = rSges(6,1) * t3542 - Icges(6,4);
t2756 = t3320 * cos(t3198) / 0.2e1 - ((2 * Icges(6,1)) - (2 * Icges(6,2)) - 0.2e1 * t3494) * sin(t3198) / 0.8e1 + (-sin(t3575) / 0.8e1 + sin(t3574) / 0.8e1) * (-Icges(6,1) + Icges(6,2) + t3494) + (-t3481 / 0.2e1 - t3480 / 0.2e1) * t3556 + ((t3632 + t2959 / 0.2e1) * rSges(6,2) + (t2954 / 0.2e1 - t2955 / 0.2e1) * rSges(6,1)) * m(6) * pkin(12) + ((-t3274 / 0.4e1 - t3275 / 0.4e1 - t3371 / 0.4e1 - t3372 / 0.4e1) * rSges(6,2) + (t3273 / 0.4e1 - t3272 / 0.4e1 - t3368 / 0.4e1 - t3369 / 0.4e1) * rSges(6,1)) * t3557 + ((t3524 / 0.4e1 + t3523 / 0.4e1 + t3522 / 0.4e1 + t3521 / 0.4e1) * rSges(6,2) + (t3519 / 0.4e1 - t3518 / 0.4e1 + t3517 / 0.4e1 - t3516 / 0.4e1) * rSges(6,1)) * t3559 + (t3597 * t3470 + t3584) * cos(t3309) / 0.16e2 - (t3596 * t3470 + t3584) * cos(t3310) / 0.16e2 - t3040 * sin(t3309) / 0.4e1 + t3041 * sin(t3310) / 0.4e1 - (cos(t3575) + cos(t3574)) * t3320 / 0.4e1;
t3465 = t2756 * qJD(1);
t3452 = (rSges(9,1) * t3484 - t3138 - t3361) * t2906;
t3450 = (m(4) * t3488 - t3188 * t3589) * t3150;
t3115 = pkin(1) * t3146;
t3024 = t3115 - t3509;
t3154 = cos(t3180);
t3449 = t3024 * t3154;
t3442 = (-rSges(4,1) * t3192 - t3489) * t3143;
t3438 = t3125 * t3130;
t3437 = t3129 * t3126;
t3434 = t3143 * (rSges(4,2) * t3550 - Icges(4,6));
t3080 = t3109 / 0.2e1;
t3081 = t3507 / 0.2e1;
t3414 = t3508 / 0.2e1 + t3077;
t3425 = (t3081 + t3080) * t3579 + t3414 * t3580;
t3317 = t3437 * t3501;
t3318 = t3438 * t3501;
t3416 = t3318 / 0.2e1 - t3317 / 0.2e1;
t3304 = t3125 * t3379;
t3305 = t3129 * t3378;
t3415 = t3147 * t3305 + t3154 * t3304;
t3396 = -0.2e1 * t3452;
t3116 = pkin(1) * t3154;
t3395 = t3081 - t3116 / 0.2e1;
t3389 = -t3512 / 0.2e1;
t3388 = -t3508 / 0.4e1;
t3387 = -t3507 / 0.4e1;
t3386 = qJD(4) * t3562;
t2894 = t3169 - t3486;
t3384 = t2894 * t3491;
t3360 = t3197 * t3150 * t3143;
t3358 = -t3449 / 0.4e1;
t3355 = -t3434 / 0.2e1;
t3354 = -t3428 / 0.2e1;
t3342 = t3355 - t3637;
t3341 = t3628 + t3355;
t2991 = -t3024 + t3512;
t3312 = t2991 * t3609;
t2939 = m(6) * t3312;
t2981 = pkin(1) * (m(4) * t3489 + t3192 * t3589) * t3143;
t3335 = t3116 - t3511;
t2992 = t3109 + t3335;
t3330 = m(6) * t3388;
t3340 = t2992 * t3330 + t2939 + t2981;
t3009 = -0.4e1 * t3173 * t3492 + 0.2e1 * t3492 + 0.2e1 * t3495;
t3336 = -t3116 - t3511;
t3334 = t3115 - 0.2e1 * t3509;
t3329 = t3485 * t3169;
t3313 = t2991 * t3387;
t3299 = t3205 * t3360;
t3287 = t3449 * t3221;
t3286 = t2861 * t3329;
t3277 = t2861 * t3316;
t3270 = t3434 / 0.2e1;
t3268 = (pkin(4) * t3548 + rSges(4,1) * t3550 - Icges(4,5)) * t3150 + t3087 * t3562 + t3088 * t3564 + (t3127 + t3128) * t3542 * t3635;
t3258 = -0.2e1 * (-t3137 * t3595 + t3076) * t3455 + t3396;
t3254 = t3631 * t3640 - t3277 - t3299;
t3237 = t3614 * pkin(4) * t3400 + t3127 * t3349 + t3128 * t3347;
t3091 = 0.2e1 * t3109;
t3235 = m(6) * (-t3334 + t3512) * t3609 + (t3091 + t3335) * t3330 + t3007 - t3299 + t3415;
t3111 = 0.4e1 * t3119;
t3105 = -t3115 / 0.2e1;
t3085 = -0.2e1 * t3101;
t2989 = 0.2e1 * t3255;
t2987 = 0.2e1 * t3322 + t3613;
t2940 = m(6) * t3313;
t2910 = -t3496 / 0.4e1;
t2895 = t3169 - 0.2e1 * t3486;
t2885 = -t3497 / 0.4e1;
t2874 = t2989 * t2907 + t3513;
t2827 = t2940 + t3287 + (pkin(1) * t3358 + t3312) * m(6);
t2787 = -0.2e1 * t3315;
t2778 = t3259 - t3586;
t2775 = (t3284 - t3626) * rSges(6,2) + (t3283 + t3282) * rSges(6,1);
t2766 = t3621 * t3044 + t3418 * t3045 - t3420 * t3046 + t3419 * t3047 + t3061 * t3152 - t3062 * t3151 + t3063 * t3145 - t3064 * t3144;
t2763 = t3270 + t3341 + t3637;
t2762 = t3270 + t3342 - t3628;
t2761 = t3268 + t3341 + t3342;
t2755 = t3009 * t3354 - t3622 + (-t2861 * t3601 + t3271 * t3631 + t3258) * m(9) + t3254 + t3416;
t2753 = (-t3360 + (t3438 / 0.4e1 - t3437 / 0.4e1) * m(6)) * t3205 + (t3329 - t3384 - t3612) * m(9) + t2787 + t3340 + t3415 + t3599 + t3600;
t2746 = t2910 + t3237;
t2745 = t3496 / 0.4e1 + t3291;
t2744 = t2910 + t3291;
t2743 = t2885 + ((-t3149 / 0.4e1 + t3148 / 0.4e1) * rSges(6,2) + (t3142 / 0.4e1 + t3141 / 0.4e1) * rSges(6,1)) * t3559 + t3237;
t2742 = t3497 / 0.4e1 + t3629;
t2741 = t2885 + t3629;
t2735 = (t3244 + 0.2e1 * t3425 + t3208) * t3633;
t1 = [-t2739 * qJD(2) - t2740 * qJD(3) + t2756 * qJD(4), -t3467 + t2761 * qJD(3) + t2741 * qJD(4) + (t3313 + (t3358 + (-t3148 / 0.2e1 - t3149 / 0.2e1) * rSges(6,2) + (t3142 / 0.2e1 - t3141 / 0.2e1) * rSges(6,1)) * pkin(1)) * t3504 + (Icges(3,5) * t3192 - t3188 * Icges(3,6) + t2939 + t3287 + t3268 + (-rSges(3,1) * t3192 + rSges(3,2) * t3188) * rSges(3,3) * m(3) + (-rSges(7,2) * t3546 + Icges(7,6)) * t3124 - (rSges(7,1) * t3546 - Icges(7,5)) * t3520 + 0.2e1 * t3355 + (-m(8) * rSges(8,3) - t3545 - t3548 - t3550) * t3169 + t2861 * t3624) * qJD(2), t2761 * qJD(2) + t2744 * qJD(4) - t3466 + (t3268 - t3434 - t3624) * qJD(3), t3465 + t2741 * qJD(2) + t2744 * qJD(3) + (t3041 * t2955 / 0.2e1 - (t3596 * m(6) - Icges(6,5)) * t2959 / 0.2e1 - t3040 * t2954 / 0.2e1 + (t3597 * m(6) - Icges(6,5)) * t3632 + (0.2e1 * pkin(12) * (-t3480 - t3481) + t3303 + t3413) * t3564 - ((-t3148 - t3149) * rSges(6,2) + (-t3141 + t3142) * rSges(6,1)) * t3559 / 0.2e1) * qJD(4); t2827 * qJD(2) + t2763 * qJD(3) + t2742 * qJD(4) + t3467, t2827 * qJD(1) + t2753 * qJD(3) + (-t3585 * t3428 + (t3024 + t3512) * t3387 + (-t3109 + t3511) * t3389 + (t3109 + t3336) * t3388) * t3504 + ((t3111 + t3009) * t3354 + 0.2e1 * t3287 + t3254 + t3340 + (-t3450 + (t3441 + t3442) * m(4)) * pkin(1) + (t2861 * t3384 - t3286 + ((t3167 * t3172 - t3107 + 0.4e1 * t3322) * t2907 - t3255 * t3462 + t3514) * t2908 - (-t3382 - t3510) * t3490 + t3396 - t3631 * t3457) * m(9) + t3595 * (-t3037 + t3120) + t3590) * qJD(2), t2763 * qJD(1) + t2753 * qJD(2) + t3468 + (t3085 + (0.16e2 * t3101 + t3112) * t3525 - 0.2e1 * t2854 + t2981 + t3235 + t3412 + ((-t3271 - 0.2e1 * t3367) * t3631 + (t3362 + t3601) * t2861 + t3329 + t2989 * t2906 + (0.4e1 * t2908 * t2987 - t2874) * t2907 + (-rSges(9,1) * t2895 + 0.2e1 * t3433) * t2903 - t3258) * m(9) + (-0.2e1 * t3289 - t3119 + t3009 / 0.2e1) * t3428 + t3622 + t3599) * qJD(3), t2742 * qJD(1) + t2736 * qJD(3) + (t2766 * t2770 * (t3169 - t3505) + t2775 * (t3512 + 0.2e1 * t3105 + t3572) + t2778 * (0.2e1 * t3530 + t3511 + 0.2e1 * t3395)) * t3386; t2762 * qJD(2) + t2745 * qJD(4) + t3466, t2762 * qJD(1) + t2755 * qJD(3) - t3468 + (t3595 * (t3085 + t3120 + 0.2e1 * t3407) + t2894 * t3402 - t3326 * t3169 + (t3102 + t3111 + (-0.8e1 * t3173 + 0.4e1) * t3492) * t3354 - t3600 + t2992 * t3126 * t3405 + t3590 - pkin(1) * t3450 + t2895 * t3328 + t2940 + t3640 + (t3091 + t3336) * t3330 + t3153 * t3304 + t3146 * t3305 - (t2989 * t3462 - t3514 + (-0.2e1 * t3167 * t3595 - 0.8e1 * t3322) * t2907) * t3498 + t2874 * t3392 - 0.4e1 * t2987 * t2907 * t3391 - t3318 / 0.4e1 + t3317 / 0.4e1 - 0.2e1 * t3277 + m(6) * (t3334 + t3512) * t3387 + t3442 * t3560 + t3235 + (t3612 - 0.4e1 * t3452 - t3286 - (-0.2e1 * t3382 - t3510) * t3490) * m(9)) * qJD(2), t2755 * qJD(2) + ((t2901 * t3540 - t3271) * m(9) + t2787 + t3416 + (-t3360 + (0.2e1 * t3174 * t3432 - t3411 * t3428) * m(5)) * t3205 + t3615) * qJD(3), t2745 * qJD(1) - t2736 * qJD(2) + (-t2766 * t2767 * t3505 + t2775 * t3572 + 0.2e1 * t2778 * (t3081 + t3530)) * t3386; t2743 * qJD(2) + t2746 * qJD(3) - t3465, t2743 * qJD(1) + ((t3080 - t3511 / 0.2e1 + t3395) * t3579 + (t3389 + t3105 + t3414) * t3580 + t2764 * t2770 * (t3506 - t3510)) * t3504 / 0.2e1 + t2735 * qJD(3), t2746 * qJD(1) + t2735 * qJD(2) + (t3208 + t3425) * qJD(3) * t3562, 0;];
Cq = t1;
