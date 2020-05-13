% Calculate inertial parameters regressor of inverse dynamics with ic joint torque vector for
% fourbarprisIC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
% 
% Output:
% tau_reg [1x(4*10)]
%   inertial parameter regressor of inverse dynamics with ic joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:59
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbarprisIC_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisIC_invdynJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisIC_invdynJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbarprisIC_invdynJ_fixb_reg2_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisIC_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisIC_invdynJ_fixb_reg2_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_regressor_ic_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:59:32
% EndTime: 2020-05-07 09:59:32
% DurationCPUTime: 0.07s
% Computational Cost: add. (45->19), mult. (71->28), div. (20->4), fcn. (70->4), ass. (0->14)
t10 = sin(qJ(1));
t12 = cos(qJ(1));
t23 = -g(1) * t10 + g(2) * t12;
t22 = (2 * qJD(2) * qJD(1)) + t23;
t21 = -g(1) * t12 - g(2) * t10;
t11 = cos(qJ(3));
t9 = sin(qJ(3));
t19 = 0.1e1 / pkin(2) / (-t10 * t11 + t9 * t12);
t18 = (-t12 * t10 - t11 * t9) / (pkin(3) + qJ(2)) / (-t11 ^ 2 + t12 ^ 2);
t17 = qJ(2) * qJDD(1);
t16 = qJDD(2) - t21;
t13 = qJD(1) ^ 2;
t1 = qJDD(1) * t18;
t2 = [0, 0, 0, 0, 0, t1, t23 * t18, t21 * t18, 0, 0, 0, 0, 0, t1, 0, 0, t16 * t18 + qJDD(1), 0, -t13 + (0.2e1 * t17 + t22) * t18, (-t13 + (t17 + t22) * t18) * qJ(2) + t16, 0, 0, 0, 0, 0, -qJDD(3) * t19, -(-g(1) * t9 + g(2) * t11) * t19, -(-g(1) * t11 - g(2) * t9) * t19, 0, 0;];
tau_reg = t2;
