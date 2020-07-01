% Calculate minimal parameter regressor of joint inertia matrix for
% palh2m1DE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% MM_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-30 17:39
% Revision: b9e8aa5c608190a7b43c48aaebfd2074f0379b0d (2020-06-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = palh2m1DE_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'palh2m1DE_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'palh2m1DE_inertiaJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-30 17:39:36
% EndTime: 2020-06-30 17:39:36
% DurationCPUTime: 0.14s
% Computational Cost: add. (47->19), mult. (109->43), div. (0->0), fcn. (113->6), ass. (0->26)
t17 = sin(qJ(2));
t19 = cos(qJ(3));
t16 = sin(qJ(3));
t20 = cos(qJ(2));
t22 = t20 * t16;
t7 = -t19 * t17 - t22;
t29 = pkin(3) * t7;
t15 = sin(qJ(4));
t13 = t19 * pkin(3) + pkin(2);
t24 = t17 * t16;
t21 = -pkin(3) * t24 + t13 * t20 + pkin(1);
t3 = pkin(4) + t21;
t28 = t15 * t3;
t27 = t16 * pkin(2);
t18 = cos(qJ(4));
t26 = t18 * t3;
t25 = t19 * pkin(2);
t23 = t17 * t20;
t14 = t20 * pkin(2) + pkin(1);
t9 = -t20 * t19 + t24;
t6 = pkin(3) * t22 + t17 * t13;
t5 = t18 * t29;
t4 = t15 * t29;
t2 = t18 * t6;
t1 = t15 * t6;
t8 = [1, 0, 0, t17 ^ 2, 0.2e1 * t23, 0, 0, 0, 0.2e1 * pkin(1) * t20, -0.2e1 * pkin(1) * t17, t7 ^ 2, 0.4e1 * (t19 ^ 2 - 0.1e1 / 0.2e1) * t23 + (0.4e1 * t20 ^ 2 - 0.2e1) * t19 * t16, 0, 0, 0, -0.2e1 * t14 * t9, 0.2e1 * t14 * t7, 0.2e1 * t21, 1, 0.2e1 * t26, -0.2e1 * t28; 0, 0, 0, 0, 0, -t17, -t20, 0, 0, 0, 0, 0, t7, t9, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t25, -0.2e1 * t27, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t9, 0, 0, 0, 0, 0, -t4, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t25, -t27, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t26, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t8;
