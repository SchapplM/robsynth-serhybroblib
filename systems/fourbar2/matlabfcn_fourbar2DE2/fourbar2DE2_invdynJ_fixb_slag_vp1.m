% Calculate vector of inverse dynamics joint torques for
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
% m [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [1x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:27
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = fourbar2DE2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(1,1),zeros(3,1),zeros(2,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar2DE2_invdynJ_fixb_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar2DE2_invdynJ_fixb_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [1 1]), ...
  'fourbar2DE2_invdynJ_fixb_slag_vp1: qJDD has to be [1x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'fourbar2DE2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2DE2_invdynJ_fixb_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2DE2_invdynJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar2DE2_invdynJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar2DE2_invdynJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:27:48
% EndTime: 2020-04-24 20:27:48
% DurationCPUTime: 0.10s
% Computational Cost: add. (59->18), mult. (160->31), div. (0->0), fcn. (64->2), ass. (0->9)
t11 = sin(qJ(1));
t12 = cos(qJ(1));
t8 = t12 * rSges(2,1) - t11 * rSges(2,2);
t6 = t11 * rSges(2,1) + t12 * rSges(2,2);
t7 = t12 * rSges(4,1) - t11 * rSges(4,2);
t5 = t11 * rSges(4,1) + t12 * rSges(4,2);
t10 = m(2) * rSges(2,2) + m(4) * rSges(4,2);
t9 = m(2) * rSges(2,1) + m(3) * pkin(2) + m(4) * rSges(4,1);
t1 = [-(-t9 * g(1) - g(2) * t10) * t11 + (t10 * g(1) - g(2) * t9) * t12 + (m(2) * (t6 ^ 2 + t8 ^ 2) + Icges(2,3) + m(3) * (t11 ^ 2 + t12 ^ 2) * pkin(2) ^ 2 + m(4) * (t5 ^ 2 + t7 ^ 2) + Icges(4,3)) * qJDD(1);];
tau = t1;
