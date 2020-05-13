% Calculate kinetic energy for
% fourbarprisOL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
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
% Datum: 2020-05-07 09:52
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbarprisOL_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'fourbarprisOL_energykin_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'fourbarprisOL_energykin_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisOL_energykin_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisOL_energykin_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisOL_energykin_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisOL_energykin_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:52:03
% EndTime: 2020-05-07 09:52:03
% DurationCPUTime: 0.03s
% Computational Cost: add. (18->17), mult. (49->35), div. (0->0), fcn. (16->4), ass. (0->12)
t90 = rSges(3,3) + qJ(2);
t87 = cos(qJ(1));
t86 = cos(qJ(3));
t85 = sin(qJ(1));
t84 = sin(qJ(3));
t83 = -t87 * rSges(2,1) + t85 * rSges(2,2);
t82 = -t86 * rSges(4,1) + t84 * rSges(4,2);
t81 = -t85 * rSges(2,1) - t87 * rSges(2,2);
t80 = -t84 * rSges(4,1) - t86 * rSges(4,2);
t79 = -qJD(2) * t85 + (-t85 * rSges(3,1) - t90 * t87) * qJD(1);
t78 = -qJD(2) * t87 + (-t87 * rSges(3,1) + t90 * t85) * qJD(1);
t1 = m(3) * (t78 ^ 2 + t79 ^ 2) / 0.2e1 + (m(4) * (t80 ^ 2 + t82 ^ 2) / 0.2e1 + Icges(4,3) / 0.2e1) * qJD(3) ^ 2 + (m(2) * (t81 ^ 2 + t83 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * qJD(1) ^ 2;
T = t1;
