% Calculate vector of inverse dynamics base forces with Newton-Euler for
% palh2m1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-03 00:53
% Revision: 702b02ffee4f47164fa6a11b998f1d39ead3f7a6 (2020-05-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = palh2m1OL_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'palh2m1OL_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'palh2m1OL_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'palh2m1OL_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'palh2m1OL_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1OL_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'palh2m1OL_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'palh2m1OL_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'palh2m1OL_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-03 00:24:34
% EndTime: 2020-05-03 00:25:48
% DurationCPUTime: 75.65s
% Computational Cost: add. (12235->1296), mult. (18961->1572), div. (0->0), fcn. (7542->10), ass. (0->715)
t3748 = cos(qJ(2));
t3743 = sin(qJ(2));
t4236 = mrSges(6,1) * pkin(6);
t3656 = -Ifges(6,5) + t4236;
t3740 = sin(qJ(5));
t3572 = t3656 * t3740;
t4233 = mrSges(6,2) * pkin(6);
t3655 = -Ifges(6,6) + t4233;
t3745 = cos(qJ(5));
t3470 = -t3655 * t3745 - t3572;
t3456 = Ifges(5,6) - t3470;
t3741 = sin(qJ(4));
t3420 = t3456 * t3741;
t4209 = mrSges(6,1) * t3740;
t3629 = pkin(4) * t4209;
t3725 = t3745 ^ 2;
t3662 = Ifges(6,4) * t3725;
t3737 = Ifges(6,1) - Ifges(6,2);
t4234 = mrSges(6,2) * pkin(4);
t4305 = (-t3737 * t3740 + t4234) * t3745;
t4342 = t4305 + Ifges(6,4) + t3629 - 0.2e1 * t3662;
t3402 = -Ifges(5,5) + t4342;
t3746 = cos(qJ(4));
t4365 = t3402 * t3746;
t3241 = t4365 + t3420;
t3747 = cos(qJ(3));
t3369 = t3402 * t3741;
t3421 = t3456 * t3746;
t3242 = t3421 - t3369;
t3742 = sin(qJ(3));
t4119 = t3742 * t3242;
t4390 = t3241 * t3747 + t4119;
t4181 = t4390 * t3743;
t4120 = t3742 * t3241;
t4386 = t3242 * t3747 - t4120;
t4395 = -t4386 * t3748 + t4181;
t4114 = t3743 * t4386;
t4394 = -t3748 * t4390 - t4114;
t4242 = pkin(6) * m(6);
t3723 = pkin(1) * t4242;
t3663 = mrSges(6,2) * t3740;
t3630 = pkin(5) * t3663;
t4206 = t3745 * mrSges(6,1);
t4244 = pkin(4) * m(6);
t3919 = mrSges(5,1) + t4206 + t4244;
t3879 = pkin(5) * t3919 - t3630;
t3872 = -t3723 + t3879;
t3739 = mrSges(5,2) - mrSges(6,3);
t4213 = t3739 * pkin(1);
t3404 = t4213 + t3872;
t3373 = t3404 * t3741;
t3756 = mrSges(4,2) * pkin(5);
t3806 = pkin(1) * m(6);
t3719 = pkin(3) * t3806;
t3804 = pkin(3) * m(5);
t3720 = pkin(1) * t3804;
t3761 = mrSges(4,1) * pkin(1);
t3970 = t3719 + t3720 + t3761;
t3933 = t3756 - t3970;
t3293 = -t3933 - t3373;
t3634 = pkin(1) * t3663;
t3899 = t3919 * pkin(1);
t3880 = t3899 - t3634;
t3661 = pkin(5) * t3739;
t4196 = pkin(5) * t4242;
t3964 = t3661 - t4196;
t3406 = t3880 - t3964;
t4161 = t3406 * t3746;
t3204 = -t4161 - t3293;
t3200 = t3204 * t3747;
t3377 = t3406 * t3741;
t3757 = mrSges(4,2) * pkin(1);
t4243 = pkin(5) * m(6);
t3718 = pkin(3) * t4243;
t3724 = pkin(5) * t3804;
t3760 = mrSges(4,1) * pkin(5);
t3971 = t3718 + t3724 + t3760;
t3935 = t3757 + t3971;
t3295 = -t3935 - t3377;
t3374 = t3404 * t3746;
t3203 = t3374 - t3295;
t4392 = -t3203 * t3742 - t3200;
t4208 = mrSges(6,2) * t3745;
t3632 = pkin(3) * t4208;
t3635 = pkin(3) * t4209;
t3754 = mrSges(5,3) * pkin(3);
t3923 = t3420 + t3632 + t3635 + t3754 - Ifges(4,5);
t3213 = t4365 + t3923;
t3206 = t3213 * t3747;
t3237 = -Ifges(4,6) - t3242;
t4391 = -t3237 * t3742 + t3206;
t4389 = t3203 * t3747 - t3204 * t3742;
t3376 = t3404 * t3743;
t3493 = -t3663 + t3919;
t3488 = pkin(2) * t3493;
t3272 = -t3488 + t3376;
t3256 = t3272 * t3746;
t3268 = t3295 * t3743;
t3664 = mrSges(6,3) + t4242;
t3626 = -mrSges(5,2) + t3664;
t4135 = t3626 * t3741;
t3477 = m(6) * pkin(3) + mrSges(4,1) + t3804 + t4135;
t3462 = pkin(2) * t3477;
t4388 = t3256 - t3268 - t3462;
t3267 = t3293 * t3743;
t4124 = t3741 * t3493;
t3443 = t4124 + mrSges(4,2);
t4220 = t3443 * pkin(2);
t4387 = -t3267 - t4220;
t4115 = t3742 * t3746;
t3201 = t3402 * t4115 + t3923 * t3742 - Ifges(3,6);
t3726 = t3746 ^ 2;
t4128 = t3726 * t3743;
t4131 = t3656 * t3745;
t4133 = t3655 * t3740;
t3471 = t4131 - t4133;
t4333 = Ifges(5,4) + t3471;
t4368 = t3664 * pkin(4);
t4372 = t4368 + t4333;
t4380 = t4372 * t4128;
t3948 = t3742 * t4380;
t4385 = t4372 * t3741;
t3458 = pkin(3) * t4124;
t4314 = Ifges(4,4) - t3458;
t3292 = t4372 - t4314;
t3487 = pkin(3) * t3493;
t4384 = t4385 + t3487 / 0.4e1;
t3473 = t3487 / 0.2e1;
t3283 = t4385 + t3473;
t4153 = t3493 * t3746;
t4295 = t3477 + t4153;
t3317 = t3742 * t4295;
t3552 = t3626 * t3746;
t3397 = -t3552 + t4124;
t3387 = mrSges(4,2) + t3397;
t3350 = t3387 * t3747;
t3211 = t3350 + t3317;
t3318 = t4295 * t3747;
t4082 = t3387 * t3742 - t3318;
t4383 = t3743 * t3211 + t3748 * t4082;
t4382 = t4372 * t3726;
t4117 = t3742 * t3726;
t4008 = t4372 * t4117;
t3534 = t4208 + t4209;
t4150 = t3534 * t3741;
t3491 = pkin(3) * t4150;
t3616 = t4233 / 0.4e1 - Ifges(6,6) / 0.4e1;
t4073 = t3616 * t3745 + t3572 / 0.4e1;
t3727 = t3747 ^ 2;
t4158 = t4342 * t3741;
t4000 = t3746 * t4158;
t4351 = t3470 * t3726;
t3893 = t4000 + t4351;
t3617 = t4233 / 0.2e1 - Ifges(6,6) / 0.2e1;
t4072 = t3617 * t3745 + t3572 / 0.2e1;
t4179 = (-t3893 - t3491 / 0.2e1 - t4072) * t3727;
t4378 = -t3491 / 0.4e1 - t4073 - t4179;
t3434 = t3470 * t3741;
t3508 = pkin(3) * t3534;
t4147 = t3534 * t3746;
t4363 = -pkin(4) * t4147 + t3434 - t3508;
t4377 = t3470 * t4115 + t4363 * t3747;
t3749 = cos(qJ(1));
t3744 = sin(qJ(1));
t4104 = t3744 * t3213;
t3755 = mrSges(5,2) * pkin(5);
t3903 = pkin(5) * t3664 - t3755;
t3394 = t3880 + t3903;
t4164 = t3394 * t3746;
t3750 = pkin(1) * mrSges(6,3);
t4235 = mrSges(5,2) * pkin(1);
t4276 = t4235 - t3750;
t3393 = t3872 + t4276;
t3358 = t3393 * t3741;
t4299 = t3358 + t3933;
t4337 = (t4164 - t4299) * t3749 - t4104;
t3166 = t4337 * t3747;
t3341 = t3369 - Ifges(4,6);
t4375 = t3534 + mrSges(4,3) + mrSges(5,3);
t3952 = t4375 * pkin(2) - Ifges(3,5);
t3198 = t3341 * t3742 - t3456 * t4115 - t3952;
t3758 = mrSges(3,2) * pkin(5);
t3762 = m(4) + m(5);
t3692 = t3762 * pkin(2);
t3640 = pkin(1) * t3692;
t3805 = pkin(2) * m(6);
t3721 = pkin(1) * t3805;
t3751 = pkin(1) * mrSges(3,1);
t3972 = t3640 + t3721 + t3751;
t3937 = -t3758 + t3972;
t4165 = t3394 * t3741;
t3290 = t4165 + t3935;
t4166 = t3393 * t3746;
t4304 = t3290 + t4166;
t4376 = -t3198 * t3744 + (t3742 * t4304 - t3937) * t3749 - t3166;
t3978 = 0.8e1 * t4333;
t3384 = -t3978 - 0.8e1 * t4368;
t3342 = t3384 * t3741;
t4327 = -0.2e1 * t4131 + 0.2e1 * t4133;
t3974 = -0.2e1 * Ifges(5,4) + t4327;
t3911 = t3974 - 0.2e1 * t4368;
t3346 = t3911 * t3741;
t3808 = m(6) * pkin(4) ^ 2;
t4127 = t3737 * t3725;
t4207 = Ifges(6,4) * t3740;
t4237 = mrSges(6,1) * pkin(4);
t4139 = (t4207 + t4237) * t3745;
t4374 = t3808 - t4127 + 0.2e1 * t4139;
t3490 = mrSges(2,2) + mrSges(3,3) + t4375;
t4099 = t3746 * t3743;
t3996 = t3742 * t4099;
t4108 = t3743 * t3742;
t4357 = -t3490 * pkin(5) - t3341 * t4108 + t3456 * t3996 + t3952 * t3743 + Ifges(2,6);
t4103 = t3744 * t3237;
t3222 = t3743 * t4103;
t3262 = t3290 * t3743;
t3360 = t3393 * t3743;
t3266 = -t3488 + t3360;
t4311 = t3266 * t3746;
t4371 = (-t3462 + t4311 + t3262) * t3749 - t3222;
t4369 = t3749 * t4304 - t4103;
t3706 = m(3) + t3762;
t4367 = t3706 * pkin(1) + mrSges(2,1);
t3386 = -t3443 + t3552;
t4343 = -mrSges(3,1) - t3805 - t3692;
t3197 = t3386 * t3742 + t3318 - t4343;
t4366 = t3197 * t3743;
t4355 = mrSges(3,2) + t3211;
t3202 = t3743 * t4355;
t3435 = t3470 * t3746;
t4362 = pkin(4) * t4150 + t3435;
t3437 = t3471 * t3746;
t3631 = pkin(4) * t3663;
t3505 = pkin(4) * t4206 + Ifges(6,3) - t3631;
t3905 = -t3505 * t3741 + t3437;
t3365 = t3394 * t3743;
t3590 = pkin(2) * t3626;
t3312 = -t3590 + t3365;
t3301 = t3312 * t3742;
t4361 = -t3487 + t3301;
t3752 = mrSges(6,3) * pkin(6);
t3708 = 0.2e1 * t3752;
t3807 = m(6) * pkin(6) ^ 2;
t4034 = Ifges(5,1) + Ifges(6,2) + t3807;
t4269 = -t4034 + Ifges(5,2) + Ifges(6,3) + t4374;
t3310 = -t4269 + t3708 + 0.2e1 * t3631;
t4106 = t3743 * t3747;
t4359 = t3213 * t4106 + t4357;
t3398 = t4153 + t4135;
t4118 = t3742 * t3398;
t3217 = t3397 * t3747 + t4118;
t4326 = t3397 * t3742 - t3398 * t3747;
t4356 = t3743 * t3217 + t4326 * t3748;
t3306 = 0.2e1 * t3310;
t4354 = t3306 * t3741;
t4156 = t3443 * t3742;
t3842 = pkin(3) ^ 2;
t3810 = m(5) * t3842;
t4322 = -Ifges(4,1) + Ifges(4,2) + t3810;
t4349 = t3306 + 0.2e1 * t4322;
t3352 = t3491 - t3470;
t3348 = t3386 * t3747;
t3353 = -t4343 - t4156;
t4267 = t3806 + t4367;
t3902 = (t3626 * t4115 + t3318 + t3353) * t3748 + t4267;
t3442 = t3477 * t3742;
t3424 = t3442 + mrSges(3,2);
t4298 = t3493 * t4115 + t3424;
t4346 = t3902 - t3743 * (-t3348 + t4298);
t3891 = t4206 / 0.4e1 + t4244 / 0.4e1 + mrSges(5,1) / 0.4e1;
t3914 = -t3630 / 0.4e1 + t3723 / 0.4e1 + t3891 * pkin(5);
t3351 = t3914 - (mrSges(5,2) / 0.4e1 - mrSges(6,3) / 0.4e1) * pkin(1);
t4345 = -t3974 + 0.2e1 * t4368;
t3976 = 0.4e1 * t4333;
t3910 = -t3976 - 0.4e1 * t4368;
t4344 = -t3978 - 0.8e1 * t4368;
t3943 = 0.4e1 * t4368 + t3976;
t3609 = -0.2e1 * t3631;
t4260 = t4269 + t3609 - 0.2e1 * t3752;
t3379 = t3406 * t3743;
t3315 = -t3590 + t3379;
t4310 = t3315 * t3742;
t4341 = -t4310 + t3487;
t4058 = t4242 / 0.2e1;
t3941 = (t4236 / 0.2e1 - Ifges(6,5) / 0.2e1) * t3745 + (t4058 + mrSges(6,3) / 0.2e1) * pkin(4) + Ifges(5,4) / 0.2e1 - t3617 * t3740;
t4340 = t3941 - Ifges(4,4) / 0.2e1;
t3589 = pkin(3) * t3626;
t4332 = t3741 * t4260;
t3250 = t4332 - t3589;
t3873 = t3723 + t3879;
t3392 = t3873 - t4276;
t3330 = t4345 * t3726;
t3890 = -t3330 + t3292;
t3866 = Ifges(3,4) + Ifges(2,5) + (t3758 + t3972) * t3743 + t3890;
t3932 = t3757 - t3971;
t3946 = -mrSges(2,1) + (-m(6) - t3706) * pkin(1);
t3395 = t3903 - t3880;
t4163 = t3395 * t3741;
t4172 = t3250 * t3746;
t4185 = (t3890 + t4172) * t3727;
t3728 = t3748 ^ 2;
t3441 = -t3458 / 0.2e1;
t3884 = t4382 + t3441 - t4340;
t3280 = t4260 * t4117;
t3803 = t3842 * m(6);
t3774 = Ifges(6,2) / 0.2e1;
t3776 = -Ifges(6,2) / 0.2e1;
t4240 = Ifges(6,1) / 0.2e1;
t3667 = t4240 + t3776;
t4259 = t4139 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1 - t3752 - t3807 / 0.2e1 + t3808 / 0.2e1 - t3667 * t3725;
t3296 = t3631 + t3774 - Ifges(6,3) / 0.2e1 - t4259;
t3536 = pkin(3) * t4135;
t4262 = t3536 - Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + t3810 / 0.2e1;
t3853 = t3296 + t4262;
t3849 = t3853 + t3803 / 0.2e1;
t3981 = t3849 * t3742 + t3280 + t4220 / 0.2e1;
t4023 = pkin(2) * t3742 / 0.2e1;
t3557 = t3589 / 0.2e1;
t4086 = t3296 * t3741 + t3557;
t3329 = t4345 * t3741;
t3252 = t3329 + t3487;
t4088 = t3252 * t3742 - t3590 / 0.2e1;
t4192 = (t4185 + (t3746 * t4088 + t3981) * t3747 + (t3493 * t4023 + t4086) * t3746 + t3477 * t4023 - Ifges(3,4) / 0.2e1 + t3884) * t3728;
t4338 = ((t3392 * t3743 - t3488) * t3742 + t3250) * t3746 + ((-t3932 + t4163) * t3743 - t3462) * t3742 + t3946 * pkin(5) + t3866 - 0.2e1 * t4185 + 0.4e1 * t4192;
t4335 = 0.4e1 * t3742;
t3236 = t3421 - t3341;
t4175 = t3236 * t3742;
t3182 = t3206 + t4175;
t3759 = mrSges(3,2) * pkin(1);
t3857 = -0.2e1 * t3536 + t4260 - t4322;
t3843 = pkin(2) ^ 2;
t3691 = t3762 * t3843;
t4291 = -t3691 - (-t3842 + t3843) * m(6);
t3846 = t3857 - Ifges(3,1) + Ifges(3,2) - t4291;
t3300 = t4260 * t3726;
t3163 = t3252 * t3746 + t3300 + t3849;
t4109 = t3743 * t3727;
t4012 = t3163 * t4109;
t3963 = 0.4e1 * t4012;
t3995 = t4260 * t4128;
t4330 = t3846 * t3743 + t3759 + t3963 - 0.2e1 * t3995;
t4329 = -t4206 / 0.2e1 - t4244 / 0.2e1;
t4278 = t4166 + t4165;
t4324 = -mrSges(3,2) + t3348 - t3317;
t4286 = -(pkin(2) * (m(6) + t3762) + mrSges(3,1)) * pkin(5) - t3759;
t4101 = t3744 * t3242;
t4320 = t3749 * t4278 + t4101;
t4102 = t3744 * t3241;
t4319 = t3749 * (t4164 - t3358) - t4102;
t4256 = 0.2e1 * pkin(1);
t4251 = 0.2e1 * t3742;
t4317 = 0.2e1 * t3743;
t4100 = t3744 * t3490;
t4316 = mrSges(1,1) - t4100;
t3461 = t3490 * t3749;
t4315 = -mrSges(1,2) - t3461;
t4113 = t3743 * t4082;
t3170 = t3211 * t3748 - t4113;
t4176 = t4326 * t3743;
t3179 = t3217 * t3748 - t4176;
t3501 = t3741 * t3747 + t4115;
t4098 = t3747 * t3746;
t4116 = t3742 * t3741;
t3502 = t4098 - t4116;
t3325 = -t3501 * t3743 + t3502 * t3748;
t4193 = qJD(5) * t3534;
t4009 = t3325 * t4193;
t4313 = qJD(1) * (qJD(3) * t3170 + qJD(4) * t3179 + t4009);
t3331 = -t3434 + t3508 / 0.2e1;
t3308 = t3331 * t4115;
t3347 = t4342 * t4117;
t4301 = -t3308 - t3347;
t4300 = 0.2e1 * t4278;
t4032 = pkin(2) * t3442;
t4297 = 0.2e1 * t4032 - 0.2e1 * Ifges(3,4);
t3446 = 0.4e1 * t3458;
t4296 = t3446 - 0.4e1 * Ifges(4,4);
t4294 = 0.2e1 * t3397;
t3515 = 0.2e1 * t3536;
t4293 = t3515 + t3810;
t4292 = 0.2e1 * t3470;
t4059 = -t4243 / 0.4e1;
t3928 = -t3634 / 0.4e1 + pkin(6) * t4059;
t3356 = t3891 * pkin(1) + t3661 / 0.4e1 + t3928;
t3936 = t3718 / 0.4e1 + t3724 / 0.4e1 - t3757 / 0.4e1 + t3760 / 0.4e1;
t4282 = -t3356 * t3741 + t3936;
t4280 = t4161 - t3373;
t4273 = Ifges(6,4) / 0.4e1 - t3662 / 0.2e1;
t4272 = Ifges(4,2) + t4293;
t3162 = t3163 * t3727;
t3231 = -t3803 + t3857;
t3216 = t3231 * t4108 / 0.2e1;
t3765 = Ifges(6,3) / 0.2e1;
t3969 = t3765 + t3776 - t3631;
t3865 = t3969 + t4259;
t3271 = t3865 * t3741;
t3559 = -t3589 / 0.2e1;
t3224 = (t3271 + t3559) * t3742;
t4225 = t3743 / 0.2e1;
t3234 = t3250 * t4225;
t3297 = t3310 * t3741;
t3249 = t3297 + t3589;
t3235 = t3249 * t3742;
t3245 = -0.2e1 * t3283 * t4108;
t3247 = t4332 + t3559;
t3382 = t4385 / 0.2e1;
t3669 = t3743 * pkin(5);
t3612 = t3669 - pkin(2);
t4134 = t3626 * t3742;
t3984 = t4134 / 0.4e1;
t3259 = (t3612 * t3984 + t3382) * t3746;
t3270 = t3296 * t3726;
t3475 = -t3487 / 0.2e1;
t3284 = -t4385 + t3475;
t3354 = 0.2e1 * t4008;
t4221 = pkin(2) * t3743;
t3614 = pkin(5) - t4221;
t4226 = t3742 / 0.4e1;
t3987 = t3614 * t4226;
t3411 = t3493 * t3987;
t4230 = t3612 / 0.4e1;
t3419 = t3493 * t4230;
t3439 = t3458 / 0.4e1;
t3474 = -t3487 / 0.4e1;
t3476 = t3488 / 0.4e1;
t4227 = -t3626 / 0.4e1;
t3492 = t3614 * t4227;
t3558 = -t3589 / 0.4e1;
t3613 = t3669 - 0.2e1 * pkin(2);
t3615 = pkin(5) - 0.2e1 * t4221;
t3625 = (m(4) / 0.4e1 + m(5) / 0.4e1) * t3843;
t3685 = t3758 / 0.4e1;
t3775 = -Ifges(6,2) / 0.4e1;
t3666 = Ifges(6,1) / 0.4e1 + t3775;
t4261 = (t4207 / 0.2e1 + t4237 / 0.2e1) * t3745 - t3631 / 0.2e1 + t3808 / 0.4e1 - t3666 * t3725;
t3279 = Ifges(6,3) / 0.4e1 - Ifges(5,1) / 0.4e1 + Ifges(5,2) / 0.4e1 + t3775 - t3807 / 0.4e1 - t3752 / 0.2e1 + t4261;
t3511 = -t3536 / 0.2e1;
t3855 = -Ifges(4,2) / 0.4e1 + Ifges(4,1) / 0.4e1 + t3511 - t3810 / 0.4e1 + t3279;
t3845 = t3855 + Ifges(3,2) / 0.4e1 - Ifges(3,1) / 0.4e1 + t3625 + (t3843 / 0.4e1 - t3842 / 0.4e1) * m(6);
t3850 = t3855 - t3803 / 0.4e1;
t3847 = t3850 + t3270;
t3862 = t3279 * t3726 - Ifges(6,3) / 0.8e1 - Ifges(5,3) / 0.8e1 + Ifges(6,2) / 0.8e1 - Ifges(5,2) / 0.8e1 - Ifges(6,1) / 0.8e1 + Ifges(5,1) / 0.8e1 - t4261;
t3851 = (t3284 * t3746 + t3847) * t3727 - Ifges(4,3) / 0.8e1 + Ifges(4,2) / 0.8e1 - Ifges(4,1) / 0.8e1 + t3862;
t3859 = t3270 - t3536 / 0.4e1 + t3279;
t3869 = (-t3742 * t3747 * t4260 - t4372) * t3726 - t4185;
t3440 = t3458 / 0.2e1;
t3894 = t3440 + t4340;
t3920 = t3440 + t4372;
t3924 = t3894 * t3742 - t4008;
t3983 = -t4116 / 0.4e1;
t3947 = t3493 * t3983;
t3985 = t4135 / 0.4e1;
t3986 = t3615 * t4226;
t3990 = -t4156 / 0.4e1;
t4021 = -t4220 / 0.2e1;
t4077 = pkin(2) * t3984 - t4385;
t4087 = -t3292 * t3742 + t3354;
t4228 = t3614 / 0.4e1;
t4229 = t3613 / 0.4e1;
t4231 = t3488 / 0.2e1;
t4268 = 0.8e1 * (t3612 * t3947 + (t3419 * t3746 + t3612 * t3985 + (-t4382 + (t3271 + t3558) * t3746 + t3439 + t3941) * t3742) * t3747 + (-(t3247 * t3746 - t3330 + t3920) * t4109 + (t3492 * t3746 + t4124 * t4228 + (-t3995 - 0.2e1 * t4384 * t4099 + (-t3536 + t4260) * t4225) * t3742) * t3747 - t4380 + (t3247 * t4225 + t3411) * t3746 + t3741 * t3614 * t3984 + t3920 * t4225) * t3748 + t3862 + (pkin(2) * t3947 + (t3354 + ((t3297 + t3557) * t3742 + t3476) * t3746 + (t3441 - t4372) * t3742 + pkin(2) * t3985) * t3747 + (t3300 + (t3329 + t3473) * t3746 + t3536 / 0.2e1 + t3296) * t3727 + t3859 + (t3474 + t4077) * t3746) * t3728 + (t3859 + (-t4385 + t3474) * t3746) * t3727 + t3259) * qJD(4) + 0.8e1 * (t3612 * t3990 + (((t3245 + t3492) * t3746 + t3216 + t3443 * t4228) * t3747 + (t3411 + t3234) * t3746 + t3477 * t3987 + (t3292 / 0.2e1 + t3869) * t3743) * t3748 + (pkin(2) * t3990 + t3162 + t3847 + ((t3235 + t3476) * t3746 + t3462 / 0.4e1 + t4087) * t3747 + (t3475 + t4077) * t3746) * t3728 + ((t3224 + t3419) * t3746 + t3477 * t4230 + t3924) * t3747 + t3851 + t3259) * qJD(3) + 0.8e1 * (t3613 * t3990 + (t3742 * t4021 + t3162 + t3845 + ((t3235 + t4231) * t3746 + t3462 / 0.2e1 + t4087) * t3747 + (t3626 * t4023 + t3284) * t3746 + t3270) * t3728 + (((t3615 * t4227 + t3245) * t3746 + t3216 + t3615 * t3443 / 0.4e1) * t3747 + (t3493 * t3986 + t3234) * t3746 + t3477 * t3986 + t3685 + (Ifges(3,4) / 0.2e1 + t3869 + t3894) * t3743) * t3748 - t4343 * t3669 / 0.4e1 + ((t3493 * t4229 + t3224) * t3746 + t3477 * t4229 + t3924) * t3747 + t3851 - t3625 - m(6) * t3843 / 0.4e1 + (t3613 * t3984 + t3382) * t3746 + Ifges(3,1) / 0.8e1 - Ifges(3,2) / 0.8e1 - Ifges(3,3) / 0.8e1) * qJD(2);
t3553 = pkin(1) * t4108 - pkin(3);
t4265 = t3553 * t3741 + t3612 * t4115;
t3165 = t3748 * t4355 + t4366;
t3194 = t3197 * t3748;
t3535 = -t3663 + t4206;
t3812 = qJD(4) ^ 2;
t3813 = qJD(3) ^ 2;
t3814 = qJD(2) ^ 2;
t3815 = qJD(1) ^ 2;
t3324 = t3501 * t3748 + t3502 * t3743;
t4013 = t3324 * t4193;
t4092 = (t4082 + t4343) * t3748 + t3202;
t4194 = qJD(4) * t4356;
t4253 = -0.2e1 * qJD(5);
t4254 = 0.2e1 * qJD(3);
t4255 = 0.2e1 * qJD(2);
t4258 = qJD(1) * t3535 * t4253 + (t4092 - t4267) * t3815 + (qJD(3) * t4383 + t4194) * t4255 + t4194 * t4254 + 0.2e1 * qJD(4) * t4013 - t3179 * qJDD(4) - t3170 * qJDD(3) - t3165 * qJDD(2) + t4356 * t3812 + t4383 * t3813 - (t3194 - t3202) * t3814;
t4257 = -0.2e1 * pkin(1);
t3826 = -0.2e1 * Ifges(6,4);
t4252 = -0.2e1 * t3443;
t4250 = 0.8e1 * t3742;
t4249 = -0.4e1 * t3744;
t4248 = -0.2e1 * t3744;
t4247 = 0.4e1 * t3744;
t4246 = 0.2e1 * t3749;
t4239 = mrSges(3,1) * pkin(5);
t4238 = mrSges(5,1) * pkin(1);
t4232 = pkin(5) * mrSges(6,3);
t4224 = -mrSges(1,3) - mrSges(2,3);
t4223 = pkin(1) * t3534;
t4222 = pkin(1) * t3743;
t3665 = m(2) + t3706;
t4211 = -m(6) - t3665;
t4203 = 0.8e1 * t4314;
t4197 = pkin(1) * t4244;
t4195 = qJD(1) * t3165;
t3303 = t3315 * t3746;
t4191 = (t4389 * t3748 + (t3303 - t4387) * t3747 - t3742 * t4388) * t3744;
t3221 = t3374 + t3377;
t4190 = ((t3221 * t3747 + t3742 * t4280) * t3748 + (-t3272 * t3741 + t3303) * t3747 + (-t3272 * t3742 - t3589) * t3746 + t4341 * t3741) * t3744;
t3537 = -pkin(1) * t3746 + pkin(5) * t3741;
t4143 = t3537 * t3749;
t3415 = t3534 * t4143;
t3464 = t3612 * t3746 + t3741 * t4222;
t3938 = -0.2e1 * t3629 + 0.4e1 * t3662 + t3826 - 0.2e1 * t4305;
t4123 = t3741 * t3612;
t3999 = t3534 * t4123;
t3538 = pkin(1) * t3741 + pkin(5) * t3746;
t4142 = t3538 * t3749;
t4002 = t3534 * t4142;
t4126 = t3740 * t3745;
t4070 = t4126 * t3826 + t4127;
t3438 = -Ifges(6,1) / 0.2e1 + t3774 + t3765 + t4070;
t4157 = t3438 * t3744;
t4005 = t3746 * t4157;
t4006 = t3741 * t4157;
t4140 = t3553 * t3746;
t4146 = t3534 * t3749;
t3150 = ((0.2e1 * t3534 * t4140 + t3938) * t3749 + ((t4002 - t4006) * t3748 - t3749 * t3999 - t3743 * t4005) * t4251 + 0.2e1 * ((t3415 + t4005) * t3748 - t3743 * t4006 + t3464 * t4146) * t3747) * qJD(5);
t4188 = t4337 * t3742;
t4184 = t4391 * t3743;
t3228 = t3236 * t3747;
t4121 = t3742 * t3213;
t4183 = (t3228 - t4121) * t3743;
t4182 = t4319 * t3742;
t3381 = t3662 + (t3667 * t3740 - t4234 / 0.2e1) * t3745 - Ifges(6,4) / 0.2e1 - t3629 / 0.2e1;
t4159 = t4342 * t3726;
t4169 = t3331 * t3746;
t4180 = (t3381 + t4159 + t4169) * t3727;
t4177 = t3201 * t3744;
t4079 = pkin(3) * t3535 + t3471 * t3741;
t3253 = t3505 * t3746 + t4079;
t4171 = t3253 * t3747;
t4170 = t3266 * t3741;
t3298 = t3310 * t3726;
t3371 = -t4127 + (0.2e1 * t4207 + t4237) * t3745 + t4240 + t3969;
t4167 = t3371 * t3741;
t3407 = t3880 + t3964;
t4160 = t3407 * t3741;
t4154 = t3493 * t3742;
t4151 = t3534 * t3614;
t4149 = t3534 * t3742;
t4148 = t3534 * t3744;
t4145 = t3535 * t3744;
t4144 = t3535 * t3749;
t4125 = t3741 * t4361;
t4122 = t3742 * t4371;
t4110 = t3743 * t3283;
t4107 = t3743 * t3744;
t4105 = t3743 * t3749;
t3579 = 0.2e1 * t3590;
t3288 = (t3579 - 0.2e1 * t3365) * t3746;
t3483 = -0.2e1 * t3488;
t3576 = 0.2e1 * t3589;
t3991 = t3743 * t4101;
t3949 = t3742 * t3991;
t3992 = t3743 * t4102;
t4097 = (((-t3749 * t4300 - 0.2e1 * t4101) * t3747 - 0.2e1 * t4182) * t3748 + ((t3288 + 0.2e1 * t4170) * t3749 + 0.2e1 * t3992) * t3747 + (((t3483 + 0.2e1 * t3360) * t3742 + t3576) * t3746 + 0.2e1 * t4125) * t3749 + 0.2e1 * t3949) * qJD(4) + t3150;
t3155 = 0.2e1 * (t3438 * t3325 * t3749 + t4342 * t3744 + (-(t3537 * t3747 + t3538 * t3742) * t3748 - t3464 * t3747 - t4140 + t3612 * t4116) * t4148) * qJD(5);
t4057 = 0.2e1 * t3746;
t3409 = t3456 * t4057;
t4096 = ((((-t3409 + 0.2e1 * t3369) * t3747 + 0.2e1 * t4120) * t3748 + 0.2e1 * t4181) * t3749 + 0.2e1 * t4190) * qJD(4) + t3155;
t3226 = t3324 * t3438 * t4253;
t3383 = 0.2e1 * t3402;
t3340 = t3383 * t3746;
t4095 = (((t3340 + 0.2e1 * t3420) * t3747 + 0.2e1 * t4119) * t3748 + 0.2e1 * t4114) * qJD(4) + t3226;
t3275 = -0.4e1 * t3310 * t4128;
t3861 = t3298 + t3865;
t4094 = 0.8e1 * (-t3803 / 0.2e1 + t3861 + (t3346 - t3487) * t3746 - t4262) * t4109 + t3275;
t3264 = 0.2e1 * t4009;
t4093 = ((t3747 * t4294 + 0.2e1 * t4118) * t3748 - 0.2e1 * t4176) * qJD(4) + t3264;
t3854 = 0.4e1 * t3536 + t4349;
t3848 = t3854 + 0.2e1 * t3803;
t4091 = t3848 * t3742 + 0.2e1 * t4220;
t4010 = t3310 * t4117;
t3277 = 0.8e1 * t4010;
t3305 = 0.4e1 * t4260;
t4090 = (-0.4e1 * t3803 - 0.4e1 * t3810 - 0.8e1 * t3536 + t3305 + 0.4e1 * Ifges(4,1) - 0.4e1 * Ifges(4,2)) * t3742 + t3277;
t3480 = -0.2e1 * t3487;
t4089 = (-t3943 * t3741 + t3480) * t3742 + t3590;
t3278 = -0.4e1 * t4010;
t4063 = t3477 * t4257;
t4085 = t3743 * t4063 + t3278;
t3285 = t3305 * t3741;
t3578 = -0.4e1 * t3589;
t4084 = t3285 + t3578;
t4083 = t3576 + t4354;
t3326 = -0.16e2 * t3948;
t4081 = pkin(1) * t4252 + t3326;
t3478 = 0.4e1 * t3487;
t3260 = -t3342 + t3478;
t4080 = t3910 * t3741 + t3480;
t4031 = pkin(2) * t4154;
t3451 = -0.2e1 * t4031;
t3577 = -0.2e1 * t3589;
t4078 = t3451 + t3577;
t4029 = pkin(2) * t4134;
t4076 = -0.2e1 * t4029 - t3342;
t4064 = 0.2e1 * t4221;
t4075 = t3626 * t4064 - t3634;
t4074 = 0.4e1 * t3470;
t3581 = t3626 * t4256;
t4065 = -0.2e1 * t4222;
t4062 = pkin(5) * t3692;
t4056 = pkin(1) * t4206;
t4049 = -0.8e1 * t4180;
t3215 = t3371 * t3746 + t4079;
t4048 = t3215 * t4246;
t4047 = t3743 * t4252;
t4044 = -0.2e1 * t4135;
t4036 = -0.2e1 * t4124;
t4035 = -0.2e1 * t3746 * (t3612 * t4149 + t3435 + t4158);
t4033 = t3626 * t4222;
t4030 = pkin(2) * t4149;
t4028 = pkin(2) * t4150;
t4026 = pkin(2) * t4156;
t4025 = t3493 * t4222;
t4024 = pkin(2) * t3534 / 0.4e1;
t4022 = t4221 / 0.2e1;
t4011 = t3292 * t4108;
t4004 = t3443 * t4108;
t3998 = t3741 * t4108;
t3993 = t3743 * t4104;
t3989 = t4351 / 0.2e1;
t3988 = t4147 / 0.4e1;
t3980 = t3381 * t3742 - t4301;
t3968 = -0.2e1 * t4033;
t3979 = t3742 * t3968 + t4083;
t3450 = 0.2e1 * t4031;
t3261 = t4299 * t3743;
t3967 = t3261 - t4220;
t3962 = t3742 * t4044;
t3961 = t3471 * t3996 - t3535 * t3614;
t3959 = t3742 * t4025;
t3958 = -t3534 * t4222 / 0.4e1;
t3957 = -t4028 / 0.4e1;
t3951 = pkin(4) * t4116 - pkin(2);
t3950 = t4116 * t4151;
t3945 = -0.4e1 * t3250 * t4108 + t3493 * t4064 + t3630 - t3723;
t3944 = t4074 * t3726 + 0.2e1 * t3491 - t4292;
t3939 = (-t3666 * t3740 + t4234 / 0.4e1) * t3745 + t3629 / 0.4e1 + t4273;
t3934 = t3756 + t3970;
t3930 = t3231 * t3742 - 0.2e1 * t3280 - t4220;
t3918 = (t4235 / 0.2e1 - t3750 / 0.2e1 + t3630 / 0.2e1 - t3723 / 0.2e1 + (-mrSges(5,1) / 0.2e1 + t4329) * pkin(5)) * t3741 - t3720 / 0.2e1 - t3719 / 0.2e1 - t3756 / 0.2e1 - t3761 / 0.2e1;
t3917 = t3351 * t3741 + t3720 / 0.4e1 + t3756 / 0.4e1 + t3761 / 0.4e1 + t3719 / 0.4e1;
t3915 = pkin(2) * t4047 + t3934;
t3907 = t3742 * t4363 - t4362 * t3747;
t3463 = -pkin(1) * t4099 + t4123;
t3906 = (-t3537 * t3742 + t3538 * t3747) * t3748 - t3463 * t3747;
t3897 = t3384 + 0.16e2 * t4382;
t3896 = t3381 * t3726 + (-t4072 * t3741 - t3508 / 0.4e1) * t3746 + t3939;
t3895 = t3384 * t3726 - t3910;
t3892 = t3746 * t3919;
t3888 = t3477 * t4064 + t3932 + 0.8e1 * t3948 - 0.4e1 * t4011;
t3448 = -0.2e1 * t3458;
t3833 = 0.2e1 * Ifges(4,4);
t3887 = t3943 * t3726 + t3448 + t3833 - t4345;
t3447 = 0.2e1 * t3458;
t3886 = t3447 + t3895;
t3885 = -t3726 * t3910 + t3448 + t3911;
t3883 = -t4238 / 0.2e1 + pkin(5) * t4058 + t3634 / 0.2e1 - t3755 / 0.2e1 + t4232 / 0.2e1 + t4329 * pkin(1);
t3878 = -t3352 / 0.2e1 - t3893;
t3877 = Ifges(6,1) + Ifges(5,3) + t3609 + t3708 + t3807 + t4374;
t3875 = t3895 + t4296;
t3871 = t3536 + t3877;
t3281 = -0.8e1 * t4332;
t3868 = ((t3281 + 0.8e1 * t3589) * t3746 + t3897 + t4203) * t3727 + t3875;
t3233 = t4084 * t3746;
t3867 = (t3233 + t3875) * t3727 + t3833 + t3885;
t3863 = Ifges(4,3) + t3877 + t4293;
t3860 = t3803 + t3863;
t3852 = t3691 + t3863 + (t3842 + t3843) * m(6) + t4286 * t3743 + Ifges(3,3);
t3844 = pkin(1) ^ 2;
t3811 = qJD(5) ^ 2;
t3605 = m(1) - t4211;
t3580 = -0.2e1 * t3590;
t3575 = 0.4e1 * t3589;
t3517 = -0.4e1 * t3536;
t3482 = 0.2e1 * t3488;
t3481 = -0.4e1 * t3487;
t3479 = 0.2e1 * t3487;
t3453 = 0.2e1 * t3462;
t3449 = -0.4e1 * t3458;
t3445 = -0.2e1 * t4025;
t3444 = t4154 * t4257;
t3430 = -0.2e1 * t4032;
t3426 = -0.2e1 * t4220;
t3403 = -t4213 + t3873;
t3400 = t3534 * (t3743 * t3951 + pkin(5));
t3372 = t3403 * t3741;
t3357 = t3392 * t3741;
t3355 = t3938 * t3742;
t3336 = -0.16e2 * t4385;
t3302 = t3312 * t3746;
t3287 = t3306 * t3726;
t3282 = 0.8e1 * t4332;
t3258 = -0.4e1 * t4110;
t3255 = t3744 * t4013;
t3254 = t3749 * t4013;
t3246 = (-t3384 - t4203) * t3742;
t3244 = (t3336 - 0.8e1 * t3487) * t3742;
t3230 = (t3282 - 0.8e1 * t3589) * t3742;
t3229 = t3237 * t3747;
t3212 = t3979 * t3746;
t3185 = t3229 + t4121;
t3178 = t3185 * t3748;
t3177 = t3229 + t3201;
t3176 = -t3201 + t3228;
t3174 = t3177 * t3748;
t3172 = t4369 * t3747;
t3168 = t3952 + t3182;
t3158 = t4324 * t3743 + t3194 + t4267;
t3152 = t3168 * t3743 + t3174;
t1 = [(-t3325 * t4144 + t4148) * t3811 - t3605 * g(1) + 0.2e1 * t3744 * t4313 + t3254 * t4254 + (t3744 * t4195 + t3254) * t4255 + (-t3325 * t4146 - t4145) * qJDD(5) + t4100 * t3815 + (-t3158 * t3744 - t3461) * qJDD(1) + t4258 * t3749; (-t3325 * t4145 - t4146) * t3811 - 0.2e1 * t3749 * t4313 + (-t3325 * t4148 + t4144) * qJDD(5) - t3461 * t3815 + (t3158 * t3749 - t4100) * qJDD(1) + (-t3749 * t4195 + t3255) * t4255 + t3255 * t4254 - t3605 * g(2) + t4258 * t3744; t4092 * qJDD(2) + t4383 * qJDD(3) + t4356 * qJDD(4) - t3605 * g(3) + t3165 * t3814 + ((((0.2e1 * mrSges(4,2) + t4294) * t3747 + 0.2e1 * t3317) * t3748 - 0.2e1 * t4113) * qJD(3) + t4093) * qJD(2) + t3170 * t3813 + t4093 * qJD(3) + t3179 * t3812 + qJD(4) * t3264 + (qJDD(5) * t3534 + t3535 * t3811) * t3324; ((((-t3437 + t4167) * t3747 + t3215 * t3742) * t3748 + t3961 + (t3215 * t3747 - t3371 * t4116) * t3743) * qJD(5) * t4248 + ((((-0.8e1 * t4000 - 0.4e1 * t3491 + t4074 - 0.8e1 * t4351) * t3727 + (t4335 * t4342 - 0.8e1 * t3308 - 0.8e1 * t3347 - 0.2e1 * t4028) * t3747 + (-0.2e1 * t4030 + 0.4e1 * t4158) * t3746 + t3944) * t3728 + 0.2e1 * (t4098 * t4151 - t3950) * t3748 + (t3944 + 0.4e1 * t4000) * t3727 + (0.4e1 * t3308 + 0.4e1 * t3347 + t3355 - 0.2e1 * t3999) * t3747 + t4035 + (t4049 + (t3352 * t4335 + 0.8e1 * t3470 * t4117 + t4000 * t4250) * t3747 + 0.4e1 * t4159 + 0.4e1 * t4169 - 0.2e1 * t4342) * t3748 * t3743) * qJD(5) + t4268) * t3749) * qJD(1) + (((t3198 - t3206) * t3748 + t3177 * t3743) * t3749 - ((-t3203 * t3748 + (-t3379 + t3579) * t3746 - t3267 + t3426) * t3742 + (-t3200 + t3937) * t3748 + t3852 + ((-t3376 + t3482) * t3746 + t3268 + t3453) * t3747 + t3479 * t3746) * t3744) * qJDD(2) + (-(t3174 + t4359) * t3744 + ((((t4213 + t3945) * t3746 + t4160 + t3888) * t3747 + ((t3661 + t4056 + t4075 + t4197 + t4238) * t3742 + t3258) * t3746 + (t3372 + t3915) * t3742 + t4330) * t3748 + (t3746 * t4083 + t3887) * t3727 + ((t3403 * t3743 - t3488) * t3742 + t3250) * t3746 + t3866 + ((-t3932 - t4160) * t3743 - t3462) * t3742 + ((t3407 * t3743 + t4089) * t3746 + (t3372 + t3934) * t3743 + t3930) * t3747 + ((t4344 * t3726 + t3233 + t3943 + t4296) * t3727 + (0.4e1 * t3280 + ((-t4344 * t3741 + t3478) * t3742 + t3580) * t3746 + t4091) * t3747 + (t3450 + t4083) * t3746 + t3887 + t4297) * t3728 + ((-t3747 * t3892 - t4115 * t4242 + t4343) * t3748 + t3946) * pkin(5)) * t3749) * qJDD(1) + ((t3178 + t4184) * t3749 + t4191) * t3813 + ((t3176 * t3748 - t4359) * t3749 + (t4192 + (t4012 + (0.2e1 * t3948 + ((-t3250 * t3742 + t4231) * t3743 - t3351) * t3746 - t4011 + t3477 * t4022 - t4282) * t3747 - t3995 / 0.2e1 + (-t4110 + (t3626 * t4022 + t3356) * t3742) * t3746 + (t3743 * t4021 + t3917) * t3742 + t3845 * t3743 + pkin(2) * t4059 - t4062 / 0.4e1 - t4239 / 0.4e1 + t3759 / 0.4e1) * t3748 + (t4086 * t3746 + t3884) * t3727 + (-t3280 / 0.2e1 + (t3284 * t3742 + t3356 * t3743 + t3590 / 0.4e1) * t3746 + t3850 * t3742 + t3917 * t3743 - t4220 / 0.4e1) * t3747 - t3941 * t3726 + ((t3351 * t3743 - t3488 / 0.4e1) * t3742 + t3279 * t3741 + t3558) * t3746 + (t4282 * t3743 - t3462 / 0.4e1) * t3742 + (t3721 / 0.4e1 + t3640 / 0.4e1 + t3685 + t3751 / 0.4e1) * t3743 + t3439 + (t4236 / 0.4e1 - Ifges(6,5) / 0.4e1) * t3745 - t3616 * t3740 + (t4242 / 0.4e1 + mrSges(6,3) / 0.4e1) * pkin(4) + pkin(1) * t4059 + Ifges(3,4) / 0.4e1 - Ifges(4,4) / 0.4e1 + Ifges(5,4) / 0.4e1 + Ifges(2,5) / 0.4e1 - t4367 * pkin(5) / 0.4e1) * t4249) * t3815 - ((t3902 - t3202) * t3744 - t4315) * g(3) + (t3152 * t3749 + t3744 * ((t3937 + t4392) * t3743 + (-t4286 + t4389) * t3748)) * t3814 + (((((t3383 * t3741 - 0.2e1 * Ifges(4,6) - t3409) * t3747 + 0.2e1 * t4121) * t3748 + 0.2e1 * t4184) * t3749 + 0.2e1 * t4191) * qJD(3) + t4096) * qJD(2) + ((t4377 * t3743 + t3907 * t3748 + t3400) * t3749 + t3471 * t3744 + (t3906 - t4265) * t4145) * t3811 + qJD(4) * t3155 + (t4395 * t3749 + t4190) * t3812 + (t4394 * t3749 - ((-t3742 * t3221 + t3747 * t4280) * t3748 + (-t3315 * t3741 - t3256) * t3747 + t4341 * t3746 + t3272 * t4116 + t3871) * t3744) * qJDD(4) + ((t3185 * t3743 - t4391 * t3748) * t3749 - (t3860 + (-t4310 + t3479) * t3746 + t4387 * t3742 - t4388 * t3747 + t4392 * t3748) * t3744) * qJDD(3) + (((t3253 * t3742 - t3747 * t3905) * t3748 + (-t3505 * t4116 + t4171) * t3743 + t3961) * t3749 + t3744 * t3352 + (-t3464 * t3742 + t3906) * t4148) * qJDD(5) + t4096 * qJD(3) - ((t3350 + t4298) * t3748 + t4211 * pkin(5) + t4366 + t4224) * g(2); (((((-0.2e1 * t3718 - 0.2e1 * t3724 - 0.2e1 * t3757 - 0.2e1 * t3760 - t4300) * t3749 + 0.2e1 * t4103) * t3747 - 0.2e1 * t4188) * t3748 + ((t4299 * t4317 + t3288 + t3426) * t3749 + 0.2e1 * t3993) * t3747 + 0.2e1 * t4122) * qJD(3) + t4097) * qJD(2) + (((t4049 * t4107 + (t4327 * t3746 + 0.2e1 * t4167) * t3749 * t3747 + t3742 * t4048 + (-t3950 / 0.2e1 + (t4169 + (t3726 - 0.1e1 / 0.2e1) * t4342) * t3743) * t4247) * t3748 - t4179 * t4247 + (t3743 * t4048 + (-t3999 / 0.2e1 + t3980) * t4247) * t3747 + (-t3371 * t3998 + t3961) * t4246) * qJD(5) + ((-0.8e1 * ((t3741 * t4024 + t3980) * t3747 + t4072 * t3726 + (t3742 * t4024 - t4158 / 0.2e1) * t3746 + t4378) * t3728 + 0.8e1 * (t3614 * t3988 - t3878 * t4108) * t3747 * t3748 + t4035) * qJD(5) + t4268) * t3744) * qJD(1) + (-t4376 * t3748 + (((-t3360 + t3482) * t3746 - t3262 + t3453) * t3749 + t3222) * t3747 + (((-t3365 + t3579) * t3742 + t3479) * t3746 + t3852 + (t3261 + t3426) * t3742) * t3749 + t3201 * t4107) * qJDD(2) + ((-t3172 + t4177 + (-t3394 * t4115 + t3742 * t4299 + t4286) * t3749) * t3748 + t4376 * t3743) * t3814 - (t3665 * pkin(5) - t3626 * t3996 + t3743 * t4343 + t4324 * t3748 - t4295 * t4106 + t4004 - t4224 + t4243) * g(1) + ((-t4103 * t3747 - t4177) * t3748 - t3993 * t3747 - t4357 * t3744 + ((((t3945 + t4276) * t3746 - t4163 - pkin(5) * t3892 + t3888) * t3747 + (t3258 + (t3755 + t4075 - t4196 - t4232 + t3899) * t3742) * t3746 + (t3357 + t3915) * t3742 - pkin(2) * t4243 - t4062 - t4239 + t4330) * t3748 + ((-t3395 * t3743 + t4089) * t3746 + (t3357 + t3934) * t3743 + t3930) * t3747 + t4338) * t3749) * t3815 + (((t3237 * t3749 + ((-t4235 / 0.4e1 + t3750 / 0.4e1 + t3914) * t3746 + (-t4197 / 0.4e1 - t4056 / 0.4e1 - t3755 / 0.4e1 + t4232 / 0.4e1 - t4238 / 0.4e1 - t3928) * t3741 + t3936) * t4249 - 0.4e1 * ((-t4153 / 0.2e1 - t3477 / 0.2e1) * pkin(2) + (t3292 - 0.2e1 * t4382 + t4172) * t3742) * t4107) * t3747 + t3201 * t3749 + (-t3759 / 0.2e1 + (t3746 * t3883 + t3918) * t3742 + (t4026 + (0.2e1 * t3283 - t4029) * t3746 + t3853 + (-t3843 / 0.2e1 + t3842 / 0.2e1) * m(6) + (-m(4) / 0.2e1 - m(5) / 0.2e1) * t3843 + t3300 + Ifges(3,1) / 0.2e1 - Ifges(3,2) / 0.2e1) * t3743 + (mrSges(3,1) / 0.2e1 + (m(6) / 0.2e1 + t3762 / 0.2e1) * pkin(2)) * pkin(5)) * t4248) * t3748 + (t3213 * t4105 + ((t3743 * t3883 + t4088) * t3746 + t3918 * t3743 + t3981) * t4248) * t3747 + t4357 * t3749 + (t3748 * t3963 + t4338) * t3744) * qJDD(1) - (-t3902 * t3749 + t4105 * t4355 - t4316) * g(3) + (((-t3744 * t3905 - t4002) * t3747 + t3742 * (t3253 * t3744 + t3415)) * t3748 + (t3253 * t4107 + t3463 * t4146) * t3747 + (t3464 * t4149 - t3352) * t3749 - t3744 * (t3505 * t3998 - t3961)) * qJDD(5) + t4097 * qJD(3) + ((-t3172 - t4188) * t3748 + ((-t3302 + t3967) * t3749 + t3993) * t3747 + t4122) * t3813 + ((-t3742 * t4369 + t3166) * t3748 - t4371 * t3747 + (t3860 + (-t3301 + t3479) * t3746 + t3967 * t3742) * t3749 + t3742 * t3993) * qJDD(3) + ((-t4320 * t3742 + t4319 * t3747) * t3748 + ((-t3312 * t3741 - t4311) * t3749 - t3991) * t3747 + (t3266 * t4116 - t4361 * t3746 + t3871) * t3749 + t3742 * t3992) * qJDD(4) + (((-t3535 * t4142 - t4362 * t3744) * t3747 - (-t3535 * t4143 - t3744 * t4363) * t3742) * t3748 + (t3463 * t4144 + t4107 * t4363) * t3747 + (t3535 * t4265 - t3471) * t3749 + (t3470 * t3996 + t3400) * t3744) * t3811 + ((-t4320 * t3747 - t4182) * t3748 + ((-t3302 + t4170) * t3749 + t3992) * t3747 + ((t3266 * t3742 + t3589) * t3746 + t4125) * t3749 + t3949) * t3812 + qJD(4) * t3150; -t4394 * t3812 + (t3182 * t3748 + t4183) * t3813 + (t3168 * t3748 + t3176 * t3743) * t3814 + (t3182 * t3743 + t3178) * qJDD(3) + t4395 * qJDD(4) + ((-0.2e1 * t4026 + (t3260 * t3746 + t3305 * t3726 + t3848) * t3727 + (t3482 * t3746 + t3453 + (0.8e1 * t4382 + (t3575 - 0.4e1 * t4332) * t3746 + 0.4e1 * Ifges(4,4) + t3449 + t3910) * t3742) * t3747 + t3846 + (0.2e1 * t4029 + t4080) * t3746 + t3287) * t3728 + (0.4e1 * (t3249 * t3746 - t3726 * t3911 - t3292) * t4109 + (t4295 * t4256 + (t3579 * t3746 + t3426 + (-0.2e1 * t3803 + t3517 + (t3481 + t3342) * t3746 + 0.4e1 * t3298 - t4349) * t3742) * t3743) * t3747 - 0.4e1 * t4380 + ((t4078 - t4354) * t3743 + t3742 * t3581) * t3746 + (0.2e1 * Ifges(3,4) - 0.2e1 * Ifges(4,4) + t3430 + t3447 - t3911) * t3743 + t3353 * t4256) * t3748 + (t3746 * t4080 + t3231 + t3287) * t3727 + (-0.4e1 * t4008 + (-t3235 + t4033) * t4057 + pkin(1) * t4047 + t3292 * t4251) * t3747 + t3300 + (t3479 - 0.2e1 * t3959 - t3346) * t3746 + t3424 * t4065 + t3708 + (t3842 + t3844) * m(6) + t3706 * t3844 + Ifges(3,1) + Ifges(2,3) + t4034 + t4070 + t4272) * qJDD(1) + t3152 * qJDD(2) + ((pkin(2) * t3535 + t4171) * t3748 + t3905 * t4106 + pkin(1) * t3535 + (-t3253 * t3743 + t3748 * t3905) * t3742) * qJDD(5) - 0.8e1 * ((t4180 + (pkin(2) * t3988 + t3742 * t3878) * t3747 + t3742 * t3957 + t3896) * t3728 + ((t4098 / 0.4e1 + t3983) * t4223 + ((t3957 - t3355 / 0.4e1 + t4301) * t3747 + t3989 - (t4030 - 0.2e1 * t4158) * t3746 / 0.4e1 - t4378) * t3743) * t3748 + t3896 * t3727 + (t3741 * t3958 + (t3989 + t4000 / 0.2e1 + t3352 / 0.4e1) * t3742) * t3747 + t3939 * t3726 + (t3742 * t3958 + t4073 * t3741 + t3508 / 0.4e1) * t3746 + t3737 * t4126 / 0.4e1 - t4273) * qJD(5) * qJD(1) + (((((t3281 + t3575) * t3746 + t3449 + t3897) * t3727 + (t3277 + ((t3336 + t3481) * t3742 + t3579) * t3746 + (t3517 + t3305) * t3742 + pkin(2) * t4036) * t3747 + (t3285 + t4078) * t3746 + pkin(2) * t3962 + t3886) * t3728 + ((pkin(1) * t4036 + t3581 * t3746 + t3326) * t3747 + t3275 + t3444 * t3746 + pkin(1) * t3962 + (t3741 * t3450 + 0.8e1 * (t3861 + t3511 + (t3346 + t3475) * t3746) * t3727 + (t3479 + t4076) * t3746 + (((t3282 + t3578) * t3742 + t3483) * t3746 + (t3446 - t3384) * t3742 + pkin(2) * t4044) * t3747 + t3515 + t3306) * t3743) * t3748 + ((t3285 + t3577) * t3746 + t3886) * t3727 + (t3278 + (t4384 * t4250 + t3445) * t3746 + t3741 * t3968 + (t3536 + t3310) * t4251) * t3747 + t3212 + 0.2e1 * t3741 * t3959 + t3885) * qJD(1) + t3226) * qJD(4) + (((((t3244 + t3579) * t3746 + t3426 + t4090) * t3747 + (t3451 + t4084) * t3746 + t3430 + t3868) * t3728 + ((((t3230 + t3483) * t3743 + t3581) * t3746 + (t3246 - 0.2e1 * t3462) * t3743 + t4081) * t3747 + ((t3478 + t4076) * t3743 + t3444) * t3746 + (t3848 + 0.2e1 * t4026) * t3743 + t3742 * t4063 + t4094) * t3748 + ((t3283 * t4250 + t3445) * t3746 + (t3803 + t3310 - Ifges(4,1) + t4272) * t4251 + t4085) * t3747 + t3212 + t4004 * t4256 + t3867) * qJD(1) + t4095) * qJD(3) - (-t4346 * t3744 + t4315) * g(1) + ((((t3340 + (0.2e1 * Ifges(5,6) - t4292) * t3741 + 0.2e1 * t3635 + 0.2e1 * t3632 + 0.2e1 * t3754 - 0.2e1 * Ifges(4,5)) * t3747 + 0.2e1 * t4175) * t3748 + 0.2e1 * t4183) * qJD(3) + (((((t3230 - 0.4e1 * t3488) * t3743 + t3581) * t3746 + (t3246 - 0.4e1 * t3462) * t3743 + t4081) * t3747 + ((-0.4e1 * t4029 + t3260) * t3743 + t3444) * t3746 + (0.2e1 * Ifges(3,1) - 0.2e1 * Ifges(3,2) + t3854 + 0.4e1 * t4026) * t3743 + t3424 * t4257 + t4291 * t4317 + t4094) * t3748 + t3353 * t4065 + (((t3244 + 0.4e1 * t3590) * t3746 - 0.4e1 * t4220 + t4090) * t3747 + (-0.4e1 * t4031 + t4084) * t3746 - 0.4e1 * t4032 + 0.4e1 * Ifges(3,4) + t3868) * t3728 + (t3450 + t3979) * t3746 + t3867 + ((t3260 * t3742 + t3445 + t3580) * t3746 + t4085 + t4091) * t3747 + t4297) * qJD(1) + t4095) * qJD(2) + ((t3951 * t3534 + t4377) * t3748 - t4223 - t3907 * t3743) * t3811 - (t4346 * t3749 + t4316) * g(2);];
tauB = t1;
