% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% fourbar1OL
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% 
% Output:
% tau_reg [4x9]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-26 17:43
% Revision: 27a48890e38af062107dd0dbc7317233bd099dca (2020-06-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbar1OL_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar1OL_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar1OL_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbar1OL_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar1OL_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1OL_invdynJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-26 17:43:19
% EndTime: 2020-06-26 17:43:19
% DurationCPUTime: 0.08s
% Computational Cost: add. (38->20), mult. (52->30), div. (0->0), fcn. (30->8), ass. (0->17)
t8 = qJD(1) + qJD(2);
t20 = qJD(1) * t8;
t19 = qJD(2) * t8;
t18 = pkin(2) * qJD(1) * qJD(2);
t11 = sin(qJ(2));
t9 = qJ(1) + qJ(2);
t5 = sin(t9);
t6 = cos(t9);
t17 = -g(1) * t5 + g(2) * t6 + t11 * t18;
t14 = cos(qJ(2));
t16 = t11 * qJDD(1) * pkin(2) - g(1) * t6 - g(2) * t5 + t14 * t18;
t15 = cos(qJ(1));
t13 = cos(qJ(3));
t12 = sin(qJ(1));
t10 = sin(qJ(3));
t7 = qJDD(1) + qJDD(2);
t1 = [qJDD(1), g(1) * t12 - g(2) * t15, g(1) * t15 + g(2) * t12, t7, (t11 * t19 + (-qJDD(1) - t7) * t14) * pkin(2) + t17, (t11 * t7 + t14 * t19) * pkin(2) + t16, 0, 0, 0; 0, 0, 0, t7, (-qJDD(1) * t14 - t11 * t20) * pkin(2) + t17, -t14 * pkin(2) * t20 + t16, 0, 0, 0; 0, 0, 0, 0, 0, 0, qJDD(3), g(1) * t10 - g(2) * t13, g(1) * t13 + g(2) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tau_reg = t1;
