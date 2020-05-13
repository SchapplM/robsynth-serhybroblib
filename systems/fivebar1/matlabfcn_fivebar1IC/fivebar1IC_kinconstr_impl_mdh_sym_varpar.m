% Implicit kinematic constraints of
% fivebar1IC
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% 
% Output:
% h [(no of constraints)x1]
%   Implicit constraint equations (e.g. closed loops)
%   In a valid robot configuration qJ this has to be zero

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:19
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function h = fivebar1IC_kinconstr_impl_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1IC_kinconstr_impl_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1IC_kinconstr_impl_mdh_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From kinconstr_impl_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:19:08
% EndTime: 2020-04-27 06:19:08
% DurationCPUTime: 0.04s
% Computational Cost: add. (15->11), mult. (8->8), div. (0->0), fcn. (8->8), ass. (0->3)
t2 = qJ(1) + qJ(2);
t1 = qJ(3) + qJ(4);
t3 = [-pkin(5) * cos(t2) + cos(qJ(1)) * pkin(2) - pkin(4) * cos(t1) - cos(qJ(3)) * pkin(3) - pkin(1); -pkin(5) * sin(t2) + sin(qJ(1)) * pkin(2) - pkin(4) * sin(t1) - sin(qJ(3)) * pkin(3); -qJ(5) - t1 + t2;];
h = t3;
