% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% fourbar2DE2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% qJDD [1x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2]';
% 
% Output:
% tau_reg [1x(1*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:27
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = fourbar2DE2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar2DE2_invdynJ_fixb_reg2_slag_vp: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar2DE2_invdynJ_fixb_reg2_slag_vp: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'fourbar2DE2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2DE2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2DE2_invdynJ_fixb_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:27:49
% EndTime: 2020-04-24 20:27:49
% DurationCPUTime: 0.05s
% Computational Cost: add. (6->3), mult. (13->6), div. (0->0), fcn. (10->2), ass. (0->5)
t4 = cos(qJ(1));
t3 = sin(qJ(1));
t2 = g(1) * t4 + g(2) * t3;
t1 = g(1) * t3 - g(2) * t4;
t5 = [0, 0, 0, 0, 0, qJDD(1), t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (qJDD(1) * pkin(2) + t1) * pkin(2), 0, 0, 0, 0, 0, qJDD(1), t1, t2, 0, 0;];
tau_reg = t5;
