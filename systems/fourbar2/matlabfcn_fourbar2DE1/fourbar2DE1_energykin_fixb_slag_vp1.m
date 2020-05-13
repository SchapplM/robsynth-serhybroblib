% Calculate kinetic energy for
% fourbar2DE1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
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
% Datum: 2020-04-24 20:22
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar2DE1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(2,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar2DE1_energykin_fixb_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar2DE1_energykin_fixb_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'fourbar2DE1_energykin_fixb_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar2DE1_energykin_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar2DE1_energykin_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar2DE1_energykin_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 20:22:53
% EndTime: 2020-04-24 20:22:53
% DurationCPUTime: 0.02s
% Computational Cost: add. (11->11), mult. (42->26), div. (0->0), fcn. (10->2), ass. (0->7)
t64 = cos(qJ(1));
t63 = sin(qJ(1));
t62 = t64 * rSges(2,1) - t63 * rSges(2,2);
t61 = t64 * rSges(4,1) - t63 * rSges(4,2);
t60 = t63 * rSges(2,1) + t64 * rSges(2,2);
t59 = t63 * rSges(4,1) + t64 * rSges(4,2);
t1 = (m(2) * (t60 ^ 2 + t62 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1 + m(4) * (t59 ^ 2 + t61 ^ 2) / 0.2e1 + Icges(4,3) / 0.2e1 + m(3) * (t63 ^ 2 + t64 ^ 2) * pkin(2) ^ 2 / 0.2e1) * qJD(1) ^ 2;
T = t1;
