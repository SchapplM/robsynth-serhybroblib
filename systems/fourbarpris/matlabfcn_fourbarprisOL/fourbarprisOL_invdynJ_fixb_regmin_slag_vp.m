% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% fourbarprisOL
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
% tau_reg [4x9]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-27 17:30
% Revision: bc59515823ab4a8d0fec19bf3bf92c32c39a66b0 (2020-06-27)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbarprisOL_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(3,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisOL_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'fourbarprisOL_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbarprisOL_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_invdynJ_fixb_regmin_slag_vp: pkin has to be [3x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-27 17:29:57
% EndTime: 2020-06-27 17:29:58
% DurationCPUTime: 0.08s
% Computational Cost: add. (16->10), mult. (30->15), div. (0->0), fcn. (16->4), ass. (0->11)
t5 = sin(qJ(1));
t7 = cos(qJ(1));
t15 = -g(1) * t5 + g(2) * t7;
t14 = (2 * qJD(2) * qJD(1)) + t15;
t12 = -g(1) * t7 - g(2) * t5;
t11 = (qJDD(1) * qJ(2));
t10 = qJDD(2) - t12;
t8 = qJD(1) ^ 2;
t6 = cos(qJ(3));
t4 = sin(qJ(3));
t1 = [qJDD(1), t15, t12, t10, (2 * t11) + t14, (t11 + t14) * qJ(2), 0, 0, 0; 0, 0, 0, qJDD(1), -t8, -qJ(2) * t8 + t10, 0, 0, 0; 0, 0, 0, 0, 0, 0, qJDD(3), -g(1) * t4 + g(2) * t6, -g(1) * t6 - g(2) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0;];
tau_reg = t1;
