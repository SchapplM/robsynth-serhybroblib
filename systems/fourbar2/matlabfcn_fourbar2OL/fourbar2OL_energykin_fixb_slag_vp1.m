% Calculate kinetic energy for
% fourbar2OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 20:32
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar2OL_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(2,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbar2OL_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbar2OL_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2OL_energykin_fixb_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2OL_energykin_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar2OL_energykin_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar2OL_energykin_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:32:17
% EndTime: 2020-04-24 20:32:17
% DurationCPUTime: 0.03s
% Computational Cost: add. (23->18), mult. (47->38), div. (0->0), fcn. (14->6), ass. (0->16)
t83 = pkin(2) * qJD(1);
t80 = cos(qJ(1));
t79 = cos(qJ(3));
t78 = sin(qJ(1));
t77 = sin(qJ(3));
t76 = qJ(1) + qJ(2);
t75 = qJD(1) + qJD(2);
t74 = cos(t76);
t73 = sin(t76);
t72 = rSges(4,1) * t79 - rSges(4,2) * t77;
t71 = rSges(4,1) * t77 + rSges(4,2) * t79;
t70 = t80 * rSges(2,1) - t78 * rSges(2,2);
t69 = t78 * rSges(2,1) + t80 * rSges(2,2);
t68 = t80 * t83 + t75 * (-t74 * rSges(3,1) + t73 * rSges(3,2));
t67 = -t78 * t83 - t75 * (-t73 * rSges(3,1) - t74 * rSges(3,2));
t1 = m(3) * (t67 ^ 2 + t68 ^ 2) / 0.2e1 + t75 ^ 2 * Icges(3,3) / 0.2e1 + (m(2) * (t69 ^ 2 + t70 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1) * qJD(1) ^ 2 + (m(4) * (t71 ^ 2 + t72 ^ 2) / 0.2e1 + Icges(4,3) / 0.2e1) * qJD(3) ^ 2;
T = t1;
